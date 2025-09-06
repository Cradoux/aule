// File exports disabled by request; keep API but no-ops when called externally.
use std::fs::File;
use std::io::Write;

use engine::world::World;

// Export a CPU reimplementation of the shader mapping for a small grid to CSV
pub fn export_raster_debug_csv(
    world: &World,
    width: u32,
    height: u32,
    path: &str,
) -> std::io::Result<()> {
    // Early return to disable writing during simulation
    if true {
        return Ok(());
    }
    let mut file = File::create(path)?;
    writeln!(file, "x,y,lon,lat,face,F,iu,iv,fu,fv,upper,upper_rule,agree,iu_plus_iv,domain_violation,tri_idx,id0,id1,id2,w0,w1,w2,depth_interp,elev_interp,elev_cpu_ref,diff")?;

    let w = width as f32;
    let h = height as f32;
    let f = world.grid.frequency;
    let (face_ids_ref, face_offs_ref) = world.grid.face_vertex_table();
    let face_ids: &[u32] = face_ids_ref;
    let face_offs: &[u32] = face_offs_ref;
    let corners = world.grid.face_corners();

    // Precompute per-face neighbor mapping across edges opposite A(0),B(1),C(2)
    // Build a map from undirected edge (lo,hi) -> (face, opp_idx)
    use std::collections::HashMap;
    let mut edge_map: HashMap<(u32, u32), (usize, usize)> = HashMap::new();
    let mut neighbor: [[u32; 3]; 20] = [[u32::MAX; 3]; 20];
    for (fid, tri) in corners.iter().enumerate().take(20) {
        for opp in 0..3usize {
            let i = (opp + 1) % 3;
            let j = (opp + 2) % 3;
            let v0 = tri[i];
            let v1 = tri[j];
            let key = if v0 < v1 { (v0, v1) } else { (v1, v0) };
            if let Some((g, gopp)) = edge_map.insert(key, (fid, opp)) {
                neighbor[fid][opp] = g as u32;
                neighbor[g][gopp] = fid as u32;
            }
        }
    }

    // Helper: row_base as in shader
    let row_base = |ii: u32| -> u32 {
        let prev = ii.saturating_sub(1);
        ii * (f + 1) - (ii * prev) / 2
    };

    for y in 0..height {
        for x in 0..width {
            let fx = x as f32 + 0.5;
            let fy = y as f32 + 0.5;
            let lon = (fx / w) * std::f32::consts::TAU - std::f32::consts::PI;
            let lat = std::f32::consts::FRAC_PI_2 - (fy / h) * std::f32::consts::PI;
            let cl = lat.cos();
            let p = [cl * lon.cos(), lat.sin(), cl * lon.sin()];

            // Face pick by max dot with per-face normal (derived from face corners)
            let mut best_f = 0usize;
            let mut best_d = -1e9f32;
            for (fi, tri) in corners.iter().enumerate().take(20) {
                let a = world.grid.pos_xyz[tri[0] as usize];
                let b = world.grid.pos_xyz[tri[1] as usize];
                let c = world.grid.pos_xyz[tri[2] as usize];
                let n = [
                    (b[1] - a[1]) * (c[2] - a[2]) - (b[2] - a[2]) * (c[1] - a[1]),
                    (b[2] - a[2]) * (c[0] - a[0]) - (b[0] - a[0]) * (c[2] - a[2]),
                    (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]),
                ];
                let d = n[0] * p[0] + n[1] * p[1] + n[2] * p[2];
                if d > best_d {
                    best_d = d;
                    best_f = fi;
                }
            }
            // Start in best face, then allow up to two one-hop neighbor rollovers
            let mut f_id = best_f as u32;
            let eps_roll = 1e-5f32;
            let mut final_face = f_id;
            // attempt up to two hops; stop early if inside
            for _ in 0..2 {
                let tri = corners[f_id as usize];
                let a = world.grid.pos_xyz[tri[0] as usize];
                let b = world.grid.pos_xyz[tri[1] as usize];
                let c = world.grid.pos_xyz[tri[2] as usize];
                // Barycentric via plane projection
                let n = normalize(cross(sub(b, a), sub(c, a)));
                let t = dot(n, sub(p, a));
                let q = sub(p, scale(n, t));
                let v0 = sub(b, a);
                let v1 = sub(c, a);
                let v2 = sub(q, a);
                let d00 = dot(v0, v0);
                let d01 = dot(v0, v1);
                let d11 = dot(v1, v1);
                let d20 = dot(v2, v0);
                let d21 = dot(v2, v1);
                let denom = (d00 * d11 - d01 * d01).max(1e-12);
                let wb_tmp = (d11 * d20 - d01 * d21) / denom;
                let wc_tmp = (d00 * d21 - d01 * d20) / denom;
                let wa_tmp = 1.0 - wb_tmp - wc_tmp;
                if wa_tmp < -eps_roll || wb_tmp < -eps_roll || wc_tmp < -eps_roll {
                    // pick the most negative component and hop across that edge
                    let mut k = 0usize; // 0->opp A, 1->opp B, 2->opp C
                    let mut vmin = wa_tmp;
                    if wb_tmp < vmin {
                        vmin = wb_tmp;
                        k = 1;
                    }
                    if wc_tmp < vmin {
                        k = 2;
                    }
                    let nf = neighbor[f_id as usize][k];
                    if nf != u32::MAX {
                        f_id = nf;
                        final_face = f_id;
                        continue;
                    }
                }
                final_face = f_id;
                break;
            }
            // Compute final barycentrics in confirmed face and clamp/renormalize
            let tri = corners[final_face as usize];
            let a = world.grid.pos_xyz[tri[0] as usize];
            let b = world.grid.pos_xyz[tri[1] as usize];
            let c = world.grid.pos_xyz[tri[2] as usize];
            let n = normalize(cross(sub(b, a), sub(c, a)));
            let t = dot(n, sub(p, a));
            let q = sub(p, scale(n, t));
            let v0 = sub(b, a);
            let v1 = sub(c, a);
            let v2 = sub(q, a);
            let d00 = dot(v0, v0);
            let d01 = dot(v0, v1);
            let d11 = dot(v1, v1);
            let d20 = dot(v2, v0);
            let d21 = dot(v2, v1);
            let denom = (d00 * d11 - d01 * d01).max(1e-12);
            let wb = (d11 * d20 - d01 * d21) / denom;
            let wc = (d00 * d21 - d01 * d20) / denom;
            let mut wa = (1.0 - wb - wc).max(0.0);
            let mut wb = wb.max(0.0);
            let wc = wc.max(0.0);
            let s = (wa + wb + wc).max(1e-9);
            wa /= s;
            wb /= s;
            let _ = wc / s;

            // Lattice mapping consistent with WGSL: u <- w_a*F (α), v <- w_b*F (β)
            let u = (wa * f as f32).clamp(0.0, f as f32);
            let v = (wb * f as f32).clamp(0.0, f as f32);
            let mut iu = u.floor() as u32;
            let mut iv = v.floor() as u32;
            let mut fu = (u - iu as f32).clamp(1e-4, 1.0 - 1e-4);
            let mut fv = (v - iv as f32).clamp(1e-4, 1.0 - 1e-4);
            if iu + iv > f - 1 {
                let s = iu + iv - (f - 1);
                if fu > fv {
                    iu -= s;
                } else {
                    iv -= s;
                }
                fu = (u - iu as f32).clamp(1e-4, 1.0 - 1e-4);
                fv = (v - iv as f32).clamp(1e-4, 1.0 - 1e-4);
            }
            let _max_u_row = (f - 1) - iv;
            let eps = 1e-6f32;
            let upper_rule = (fu + fv) >= (1.0 - eps);
            let upper = upper_rule; // authoritative
                                    // Tri index
                                    // Canonical tri index per senior dev: tri_idx = face*F*F + v*(2F - v) + 2*u + (upper?1:0)
            let v = iv as u32;
            let u = iu as u32;
            let tri_local = v * (2 * f - v) + 2 * u + if upper { 1 } else { 0 };
            let tri_idx = f_id * f * f + tri_local;
            // Vertex ids for current cell
            let base = face_offs[f_id as usize];
            let id_v0 = face_ids[(base + row_base(iu) + iv) as usize];
            let id_v1 = face_ids[(base + row_base(iu + 1) + iv) as usize];
            let id_v2 = face_ids[(base + row_base(iu) + (iv + 1)) as usize];
            let id_v3 = face_ids[(base + row_base(iu + 1) + (iv + 1)) as usize];
            let (id0, id1, id2, w0, w1, w2) = if upper {
                // Upper: with u=α, v=β → V1 gets fu, V2 gets fv
                (id_v1, id_v3, id_v2, fu, 1.0 - fu - fv, fv)
            } else {
                // Lower: with u=α, v=β → V1 gets fu, V2 gets fv
                (id_v0, id_v1, id_v2, 1.0 - fu - fv, fu, fv)
            };
            // Final clamp+normalize to match shader
            let (mut ww0, mut ww1, mut ww2) = (w0.max(0.0), w1.max(0.0), w2.max(0.0));
            let s_w = (ww0 + ww1 + ww2).max(1e-9);
            ww0 /= s_w;
            ww1 /= s_w;
            ww2 /= s_w;
            let depth = ww0 * world.depth_m[id0 as usize]
                + ww1 * world.depth_m[id1 as usize]
                + ww2 * world.depth_m[id2 as usize];
            let elev = world.sea.eta_m - depth;

            // CPU reference color elevation (nearest vertex)
            let elev_cpu = elev; // placeholder: CPU painter uses nearest; here we keep interp for diff=0
            let diff = elev - elev_cpu;

            let iu_plus_iv = (iu + iv) as u32;
            let domain_violation = if iu_plus_iv > (f - 1) { 1 } else { 0 };
            let agree = if upper == upper_rule { 1 } else { 0 };
            writeln!(
                file,
                "{},{},{:.6},{:.6},{},{},{},{},{},{:.6},{:.6},{},{},{},{},{},{},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                x, y, lon, lat, f_id, f, iu, iv, fu, fv, upper as u32, upper_rule as u32, agree, iu_plus_iv, domain_violation, tri_idx, id0, id1, id2, ww0, ww1, ww2, depth, elev, elev_cpu, diff
            )?;
        }
    }
    Ok(())
}

fn sub(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn dot(a: [f32; 3], b: [f32; 3]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn cross(a: [f32; 3], b: [f32; 3]) -> [f32; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn length(a: [f32; 3]) -> f32 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn normalize(a: [f32; 3]) -> [f32; 3] {
    let l = length(a).max(1e-9);
    [a[0] / l, a[1] / l, a[2] / l]
}
fn scale(a: [f32; 3], s: f32) -> [f32; 3] {
    [a[0] * s, a[1] * s, a[2] * s]
}
