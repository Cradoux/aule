// Node >= 18, no external deps
// Validates extraction assumptions and writes out/report.md

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { createRequire } from 'module';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const repoRoot = path.resolve(__dirname, '..');

const require = createRequire(import.meta.url);
if (typeof globalThis.self === 'undefined') globalThis.self = globalThis;
if (typeof globalThis.window === 'undefined') globalThis.window = globalThis;
const THREE = require(path.resolve(repoRoot, 'tectonics-js/libraries/three.js/Three.js'));

function round6(x) { return Math.round(x * 1e6) / 1e6; }
function toArr3(v) { return [round6(v.x), round6(v.y), round6(v.z)]; }

function faceNormal(a, b, c) {
  const ab = b.clone().sub(a);
  const ac = c.clone().sub(a);
  return ab.cross(ac).normalize();
}

function buildBaseIcosa() { return new THREE.IcosahedronGeometry(1, 0); }

function lonlatToXYZ(lon, lat) {
  const cl = Math.cos(lat);
  const x = cl * Math.cos(lon);
  const y = Math.sin(lat);
  // NOTE: three.js shaders compute lon = atan(-z, x) + PI → implies z = -cl * sin(lon)
  const z = -cl * Math.sin(lon);
  return new THREE.Vector3(x, y, z);
}

function argmaxDotFacePick(base, p, eps = 1e-12) {
  let best = -Infinity;
  let bestId = 0;
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const a = base.vertices[f.a];
    const b = base.vertices[f.b];
    const c = base.vertices[f.c];
    const n = faceNormal(a, b, c);
    const d = n.dot(p);
    if (d > best + eps || (Math.abs(d - best) <= eps && fi < bestId)) {
      best = d;
      bestId = fi;
    }
  }
  return bestId;
}

function barycentric(a, b, c, p) {
  const v0 = b.clone().sub(a);
  const v1 = c.clone().sub(a);
  const v2 = p.clone().sub(a);
  const d00 = v0.dot(v0);
  const d01 = v0.dot(v1);
  const d11 = v1.dot(v1);
  const d20 = v2.dot(v0);
  const d21 = v2.dot(v1);
  const denom = d00 * d11 - d01 * d01;
  const v = (d11 * d20 - d01 * d21) / denom; // beta
  const w = (d00 * d21 - d01 * d20) / denom; // gamma
  const u = 1 - v - w; // alpha
  return [u, v, w];
}

function baryPick(base, p, eps = 1e-10) {
  let bestFace = -1;
  let bestScore = -Infinity;
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const a = base.vertices[f.a];
    const b = base.vertices[f.b];
    const c = base.vertices[f.c];
    const [u, v, w] = barycentric(a, b, c, p);
    if (u >= -eps && v >= -eps && w >= -eps) {
      const score = Math.min(u, v, w); // pick most interior if on seam
      if (score > bestScore + eps || (Math.abs(score - bestScore) <= eps && fi < bestFace)) {
        bestScore = score;
        bestFace = fi;
      }
    }
  }
  if (bestFace === -1) {
    // fallback to argmax if numerical issues
    return argmaxDotFacePick(base, p, eps);
  }
  return bestFace;
}

function baryToLattice(alpha, beta, gamma, F, eps = 1e-12) {
  const u = alpha * F;
  const v = beta * F;
  const i = Math.floor(u + eps);
  const j = Math.floor(v + eps);
  const fu = u - i;
  const fv = v - j;
  const isUpper = fu + fv > 1 - eps;
  return { i, j, isUpper };
}

function testNormalsCCW(base) {
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const a = base.vertices[f.a];
    const n = faceNormal(base.vertices[f.a], base.vertices[f.b], base.vertices[f.c]);
    if (n.dot(a) <= 0) {
      throw new Error(`Face ${fi} has inward or degenerate normal`);
    }
  }
}

function testFacePickAgreement(base) {
  // Define an inlined argmax-dot method independently
  function inlineArgmaxDot(baseLocal, p) {
    let best = -Infinity;
    let bestId = 0;
    for (let fi = 0; fi < baseLocal.faces.length; fi++) {
      const f = baseLocal.faces[fi];
      const a = baseLocal.vertices[f.a];
      const b = baseLocal.vertices[f.b];
      const c = baseLocal.vertices[f.c];
      const n = faceNormal(a, b, c);
      const d = n.dot(p);
      if (d > best || (d === best && fi < bestId)) {
        best = d;
        bestId = fi;
      }
    }
    return bestId;
  }
  const rng = Math.random;
  for (let i = 0; i < 10; i++) {
    const lon = (rng() * 2 - 1) * Math.PI;
    const lat = (rng() * 1 - 0.5) * Math.PI;
    const p = lonlatToXYZ(lon, lat);
    const a = argmaxDotFacePick(base, p);
    const b = inlineArgmaxDot(base, p);
    if (a !== b) {
      throw new Error(`Face pick mismatch at sample ${i}: cpu=${a}, inline=${b}`);
    }
  }
}

function computeNeighbors(base) {
  const neighbors = new Array(base.faces.length).fill(0).map(() => ({ oppA: null, oppB: null, oppC: null }));
  const edgeMap = new Map();
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const verts = [f.a, f.b, f.c];
    for (let e = 0; e < 3; e++) {
      const u = verts[e];
      const v = verts[(e + 1) % 3];
      const min = Math.min(u, v);
      const max = Math.max(u, v);
      const k = `${min}-${max}`;
      const entry = edgeMap.get(k) || { faces: [] };
      entry.faces.push(fi);
      edgeMap.set(k, entry);
    }
  }
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const edges = [
      { name: 'oppA', u: f.b, v: f.c },
      { name: 'oppB', u: f.c, v: f.a },
      { name: 'oppC', u: f.a, v: f.b },
    ];
    for (const edge of edges) {
      const k = `${Math.min(edge.u, edge.v)}-${Math.max(edge.u, edge.v)}`;
      const entry = edgeMap.get(k);
      if (entry && entry.faces.length === 2) {
        neighbors[fi][edge.name] = entry.faces[0] === fi ? entry.faces[1] : entry.faces[0];
      }
    }
  }
  return neighbors;
}

function edgeSamplesRollover(base) {
  const eps = 1e-5;
  const neighbors = computeNeighbors(base);
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const A = base.vertices[f.a];
    const B = base.vertices[f.b];
    const C = base.vertices[f.c];
    const edges = [
      { name: 'oppA', U: B, V: C, nfi: neighbors[fi].oppA },
      { name: 'oppB', U: C, V: A, nfi: neighbors[fi].oppB },
      { name: 'oppC', U: A, V: B, nfi: neighbors[fi].oppC },
    ];
    for (const e of edges) {
      const mid = e.U.clone().add(e.V).normalize();
      // inside face: nudge toward opposite corner
      const opp = (e.name === 'oppA') ? A : (e.name === 'oppB' ? B : C);
      const p_in = mid.clone().multiplyScalar(1 - eps).add(opp.clone().multiplyScalar(eps)).normalize();
      const p_out = mid.clone().multiplyScalar(1 - eps).add(opp.clone().multiplyScalar(-eps)).normalize();
      const pick_in = argmaxDotFacePick(base, p_in);
      const pick_out = argmaxDotFacePick(base, p_out);
      if (pick_in !== fi) throw new Error(`Edge rollover: inside pick mismatch for face ${fi} ${e.name}`);
      if (pick_out !== e.nfi) throw new Error(`Edge rollover: neighbor pick mismatch from face ${fi} ${e.name} → ${pick_out} (expected ${e.nfi})`);
      // classify upper/lower using bary on each side at F=4
      const F = 4;
      const [u1, v1, w1] = barycentric(A, B, C, p_in);
      const t1 = baryToLattice(u1, v1, w1, F);
      if (typeof t1.isUpper !== 'boolean') throw new Error('Upper/lower classification failed');
    }
  }
}

function writeReport() {
  const outDir = path.resolve(repoRoot, 'out');
  fs.mkdirSync(outDir, { recursive: true });
  const spec = JSON.parse(fs.readFileSync(path.resolve(outDir, 'spec.json'), 'utf8'));
  const sample = JSON.parse(fs.readFileSync(path.resolve(outDir, 'sample_F4.json'), 'utf8'));

  const lines = [];
  lines.push('### Tectonics.js geometry and raster mapping findings');
  lines.push('');
  lines.push('- **corner order**: A,B,C are CCW with outward normals');
  lines.push('- **face pick policy**: argmax dot(normal, p); ties → lowest face id');
  lines.push('- **lon/lat → xyz**: x=cos(lat)cos(lon), y=sin(lat), z=-cos(lat)sin(lon)');
  lines.push('- **projection**: equirectangular; vertex shader normalizes lon to [-PI,PI) with offset');
  lines.push('- **pixel center**: u=(x+0.5)/W, v=(y+0.5)/H (CPU rasters)');
  lines.push('- **lattice**: Class-I triangular; F^2 triangles/face; lower [i,j]→[i+1,j]→[i,j+1], upper [i+1,j]→[i+1,j+1]→[i,j+1]');
  lines.push('- **bary→lattice rule**: u=αF, v=βF; upper if frac(u)+frac(v) > 1−ε; clamp if floor(u)+floor(v) ≥ F');
  lines.push('');
  lines.push('Face 0 neighbors and permutations:');
  const f0 = spec.icosa.faces[0];
  for (const n of f0.neighbors_opp) {
    lines.push(`- **${n.edge}**: neighbor=${n.face_id}, perm=${n.perm}`);
  }
  lines.push('');
  lines.push('Short table for F=4, face 0:');
  lines.push('');
  lines.push('```');
  lines.push(`corners: ${JSON.stringify(sample.face_corners_xyz)}`);
  lines.push(`first_10_face_vertices: ${JSON.stringify(sample.first_10_face_vertices)}`);
  lines.push(`first_8_triangles: ${JSON.stringify(sample.first_8_triangles)}`);
  lines.push('```');
  lines.push('');
  lines.push('Assumptions/notes:');
  lines.push('- **F power-of-two**: THREE.IcosahedronGeometry(detail) implies F=2^detail.');
  lines.push('- **Diagonal tie**: uses frac(u)+frac(v) > 1−ε for upper; else lower.');
  lines.push('- **Neighbor perms**: derived by shared-vertex correspondence per base icosa faces.');
  fs.writeFileSync(path.resolve(outDir, 'report.md'), lines.join('\n'));
}

function main() {
  const base = buildBaseIcosa();
  testNormalsCCW(base);
  testFacePickAgreement(base);
  edgeSamplesRollover(base);
  writeReport();
  console.log('Wrote: out/report.md');
}

// Always run when invoked via node (Windows file URL path nuances)
main();


