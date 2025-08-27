// Node >= 18, no external deps
// Extracts canonical geometry and mapping spec from THREE.IcosahedronGeometry
// and emits out/spec.json and out/sample_F4.json

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import { createRequire } from 'module';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const repoRoot = path.resolve(__dirname, '..');

const require = createRequire(import.meta.url);

// Minimal browser shims for three.js bundling quirks
if (typeof globalThis.self === 'undefined') globalThis.self = globalThis;
if (typeof globalThis.window === 'undefined') globalThis.window = globalThis;

const THREE = require(path.resolve(repoRoot, 'tectonics-js/libraries/three.js/Three.js'));

function toArr(v) {
  return [round6(v.x), round6(v.y), round6(v.z)];
}

function round6(x) {
  return Math.round((x + 0) * 1e6) / 1e6;
}

function vec(x, y, z) {
  return new THREE.Vector3(x, y, z);
}

function normalize(v) {
  const out = v.clone();
  out.normalize();
  return out;
}

function faceNormal(a, b, c) {
  const ab = b.clone().sub(a);
  const ac = c.clone().sub(a);
  return ab.cross(ac).normalize();
}

function buildIcosaGeometry(F) {
  if (F <= 0) {
    throw new Error(`F must be positive; got ${F}`);
  }
  const detail = Math.log2(F);
  if (Math.abs(detail - Math.round(detail)) > 1e-9) {
    throw new Error(`F must be a power of 2 to match THREE.IcosahedronGeometry(detail). Got F=${F}`);
  }
  return new THREE.IcosahedronGeometry(1, Math.round(detail));
}

function buildBaseIcosa() {
  return new THREE.IcosahedronGeometry(1, 0);
}

function computeBaseNeighbors(baseGeom) {
  // Return neighbors per face across each opposite corner edge (opp A=BC, opp B=CA, opp C=AB)
  // Using face indices from base geometry (detail=0)
  const neighbors = new Array(baseGeom.faces.length).fill(0).map(() => ({ oppA: null, oppB: null, oppC: null }));
  const edgeMap = new Map(); // key: "min-max" -> {faces: [fi1, fi2], verts:[min,max]}
  for (let fi = 0; fi < baseGeom.faces.length; fi++) {
    const f = baseGeom.faces[fi];
    const verts = [f.a, f.b, f.c];
    for (let e = 0; e < 3; e++) {
      const u = verts[e];
      const v = verts[(e + 1) % 3];
      const min = Math.min(u, v);
      const max = Math.max(u, v);
      const k = `${min}-${max}`;
      const entry = edgeMap.get(k) || { faces: [], verts: [min, max] };
      entry.faces.push(fi);
      edgeMap.set(k, entry);
    }
  }
  // For each face, identify across which edge we find a neighbor
  for (let fi = 0; fi < baseGeom.faces.length; fi++) {
    const f = baseGeom.faces[fi];
    const A = f.a, B = f.b, C = f.c;
    const edges = [
      { name: 'oppA', u: B, v: C },
      { name: 'oppB', u: C, v: A },
      { name: 'oppC', u: A, v: B },
    ];
    for (const edge of edges) {
      const min = Math.min(edge.u, edge.v);
      const max = Math.max(edge.u, edge.v);
      const k = `${min}-${max}`;
      const entry = edgeMap.get(k);
      if (!entry || entry.faces.length < 2) continue;
      const nfi = entry.faces[0] === fi ? entry.faces[1] : entry.faces[0];
      neighbors[fi][edge.name] = nfi;
    }
  }
  return neighbors;
}

function neighborPermutation(ourFace, neighborFace, across) {
  // across: 'oppA' means across BC; 'oppB' across CA; 'oppC' across AB
  // Define mapping (A,B,C)->(?, ?, ?) where RHS are neighbor letters 'A','B','C'
  const ourVerts = [ourFace.a, ourFace.b, ourFace.c];
  const neighVerts = [neighborFace.a, neighborFace.b, neighborFace.c];
  const ourLetters = ['A', 'B', 'C'];
  const neighLetters = ['A', 'B', 'C'];

  const ourLetterOf = new Map([
    [ourVerts[0], 'A'],
    [ourVerts[1], 'B'],
    [ourVerts[2], 'C'],
  ]);
  const neighLetterOf = new Map([
    [neighVerts[0], 'A'],
    [neighVerts[1], 'B'],
    [neighVerts[2], 'C'],
  ]);

  // Determine shared vertices and map them
  const map = { A: null, B: null, C: null };
  const usedNeigh = new Set();
  for (const v of ourVerts) {
    if (neighVerts.includes(v)) {
      const ourL = ourLetterOf.get(v);
      const neighL = neighLetterOf.get(v);
      map[ourL] = neighL;
      usedNeigh.add(neighL);
    }
  }
  // The remaining our letter maps to the neighbor letter not used yet
  const remainingOur = ourLetters.find(L => map[L] === null);
  const remainingNeigh = neighLetters.find(L => !usedNeigh.has(L));
  if (remainingOur && remainingNeigh) {
    map[remainingOur] = remainingNeigh;
  }

  return `(A,B,C)->(${map.A}',${map.B}',${map.C}')`;
}

function lonlatToXYZ(lon, lat) {
  // Using conventional spherical to Cartesian mapping as reported in spec
  const cl = Math.cos(lat);
  const x = cl * Math.cos(lon);
  const y = Math.sin(lat);
  const z = cl * Math.sin(lon);
  return [x, y, z];
}

function pickFaceByDot(baseGeom, p, eps = 1e-12) {
  // Argmax of dot(normal, p). For base faces (detail=0)
  let best = -Infinity;
  let bestId = 0;
  for (let fi = 0; fi < baseGeom.faces.length; fi++) {
    const f = baseGeom.faces[fi];
    const a = baseGeom.vertices[f.a];
    const b = baseGeom.vertices[f.b];
    const c = baseGeom.vertices[f.c];
    const n = faceNormal(a, b, c);
    const d = n.x * p.x + n.y * p.y + n.z * p.z;
    if (d > best + eps || (Math.abs(d - best) <= eps && fi < bestId)) {
      best = d;
      bestId = fi;
    }
  }
  return bestId;
}

function barycentric(a, b, c, p) {
  // Compute barycentric coords of p in triangle abc on the plane (not reprojected)
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

function baryToLattice(alpha, beta, gamma, F, eps = 1e-12) {
  // Class-I triangular lattice on [A,B,C]
  const u = alpha * F;
  const v = beta * F;
  let i = Math.floor(u + eps);
  let j = Math.floor(v + eps);
  const fu = u - i;
  const fv = v - j;
  const isUpper = fu + fv > 1 - eps;
  // Reflect across diagonal if in upper triangle
  if (i + j >= F) {
    // clamp into valid range
    i = Math.max(0, Math.min(F, i));
    j = Math.max(0, Math.min(F - i, j));
  }
  return { i, j, isUpper };
}

function buildFaceVertexTable(geomF, baseGeom, baseFaceId, F, eps = 1e-6) {
  // Reconstruct v[i][j] for the given base face using the same lerp scheme as THREE
  const f0 = baseGeom.faces[baseFaceId];
  const A = baseGeom.vertices[f0.a];
  const B = baseGeom.vertices[f0.b];
  const C = baseGeom.vertices[f0.c];

  const table = [];
  for (let i = 0; i <= F; i++) {
    const row = [];
    const aj = A.clone().lerp(C, i / F);
    const bj = B.clone().lerp(C, i / F);
    const rows = F - i;
    for (let j = 0; j <= rows; j++) {
      const p = aj.clone().lerp(bj, rows === 0 ? 0 : j / rows);
      p.normalize();
      // Find nearest vertex id in geomF
      const vid = findNearestVertexId(geomF, p, eps);
      row.push(vid);
    }
    table.push(row);
  }
  return table;
}

function findNearestVertexId(geom, p, eps = 1e-6) {
  let best = Infinity;
  let bestId = -1;
  for (let i = 0; i < geom.vertices.length; i++) {
    const v = geom.vertices[i];
    const dx = v.x - p.x;
    const dy = v.y - p.y;
    const dz = v.z - p.z;
    const d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < best) {
      best = d2;
      bestId = i;
    }
  }
  if (best > eps * eps) {
    // Still return bestId, but flagging could be helpful during debugging
  }
  return bestId;
}

function buildTrianglesForFace(table, baseFaceId, F) {
  // Emit triangle indices per THREE's subdivision order for a single base face
  // Using the per-(i,j) vertex ids in 'table'
  const tris = [];
  for (let i = 0; i < F; i++) {
    const rows = F - i;
    for (let k = 0; k < rows; k++) {
      // lower triangle (even j): (v[i][k+1], v[i+1][k], v[i][k])
      tris.push({
        kind: 'lower',
        idx: [table[i][k + 1], table[i + 1][k], table[i][k]],
      });
      // upper triangle (odd j): (v[i][k+1], v[i+1][k+1], v[i+1][k])
      if (k + 1 <= rows - 1) {
        tris.push({
          kind: 'upper',
          idx: [table[i][k + 1], table[i + 1][k + 1], table[i + 1][k]],
        });
      }
    }
  }
  return tris;
}

function main() {
  fs.mkdirSync(path.resolve(repoRoot, 'out'), { recursive: true });

  const base = buildBaseIcosa();
  const neighbors = computeBaseNeighbors(base);

  // Assemble faces spec
  const facesSpec = [];
  for (let fi = 0; fi < base.faces.length; fi++) {
    const f = base.faces[fi];
    const A = base.vertices[f.a];
    const B = base.vertices[f.b];
    const C = base.vertices[f.c];
    const n = faceNormal(A, B, C);
    const neighIds = [neighbors[fi].oppA, neighbors[fi].oppB, neighbors[fi].oppC];
    const neighPerms = [
      neighborPermutation(f, base.faces[neighbors[fi].oppA], 'oppA'),
      neighborPermutation(f, base.faces[neighbors[fi].oppB], 'oppB'),
      neighborPermutation(f, base.faces[neighbors[fi].oppC], 'oppC'),
    ];
    facesSpec.push({
      id: fi,
      corner_order: ['A', 'B', 'C'],
      corners_xyz_ccw_outward: true,
      corners: [toArr(A), toArr(B), toArr(C)],
      normal: toArr(n),
      neighbors_opp: [
        { edge: 'BC (opp A)', face_id: neighIds[0], perm: neighPerms[0] },
        { edge: 'CA (opp B)', face_id: neighIds[1], perm: neighPerms[1] },
        { edge: 'AB (opp C)', face_id: neighIds[2], perm: neighPerms[2] },
      ],
    });
  }

  const spec = {
    icosa: {
      faces: facesSpec,
      tie_breaks: {
        face_pick: 'argmax dot(n, p) | tie→lowest face id',
        upper_tri_rule: 'fu+fv>(1-eps) picks upper within cell; otherwise lower',
      },
    },
    mapping: {
      lonlat_to_xyz: 'x=cos(lat)*cos(lon), y=sin(lat), z=-cos(lat)*sin(lon)',
      lon_range: '[-π, π]',
      lat_range: '[-π/2, π/2]',
      pixel_center: 'u=(x+0.5)/W, v=(y+0.5)/H',
      projection: 'equirectangular',
      face_pick_policy: 'argmax dot(normal, p) with epsilon E=1e-12',
    },
    lattice: {
      class: 'Class-I triangular',
      indexing: {
        row_linearization: 'row_base(i)=i*(F+1) - (i*(i-1))/2',
        idx: 'Per-row packed; see triangles order in THREE.PolyhedronGeometry',
      },
      bary_to_lattice: {
        alpha_beta: 'u=α*F, v=β*F',
        reflect: 'if floor(u)+floor(v) >= F then clamp into domain',
        which_upper: 'isUpper when frac(u)+frac(v) > 1-eps',
        epsilons: { bary: '1e-12', reflect: '1e-12' },
      },
      triangles: {
        per_face: 'F^2',
        lower_order: '[i,j]→[i+1,j]→[i,j+1]',
        upper_order: '[i+1,j]→[i+1,j+1]→[i,j+1]',
        winding_ccw_outward: true,
      },
    },
  };

  const outSpecPath = path.resolve(repoRoot, 'out/spec.json');
  fs.writeFileSync(outSpecPath, JSON.stringify(spec, null, 2));

  // Build ground truth for F=4, face 0
  const F = 4;
  const geomF = buildIcosaGeometry(F);
  const fvt = buildFaceVertexTable(geomF, base, 0, F);
  const tris = buildTrianglesForFace(fvt, 0, F);

  const sample = {
    F,
    face_id: 0,
    face_corners_xyz: spec.icosa.faces[0].corners,
    first_10_face_vertices: [],
    first_8_triangles: [],
    neighbors_opp: spec.icosa.faces[0].neighbors_opp,
  };
  // Dump first 10 (i,j) entries linearized
  outer: for (let i = 0; i < fvt.length; i++) {
    for (let j = 0; j < fvt[i].length; j++) {
      sample.first_10_face_vertices.push({ i, j, id: fvt[i][j] });
      if (sample.first_10_face_vertices.length >= 10) break outer;
    }
  }
  for (let t = 0; t < Math.min(8, tris.length); t++) {
    sample.first_8_triangles.push({ kind: tris[t].kind, idx: tris[t].idx });
  }
  const outSamplePath = path.resolve(repoRoot, 'out/sample_F4.json');
  fs.writeFileSync(outSamplePath, JSON.stringify(sample, null, 2));

  console.log(`Wrote: ${path.relative(repoRoot, outSpecPath)} and ${path.relative(repoRoot, outSamplePath)}`);
}

// Always run when invoked via node
main();


