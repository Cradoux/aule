### Tectonics.js geometry and raster mapping findings

- **corner order**: A,B,C are CCW with outward normals
- **face pick policy**: argmax dot(normal, p); ties → lowest face id
- **lon/lat → xyz**: x=cos(lat)cos(lon), y=sin(lat), z=-cos(lat)sin(lon)
- **projection**: equirectangular; vertex shader normalizes lon to [-PI,PI) with offset
- **pixel center**: u=(x+0.5)/W, v=(y+0.5)/H (CPU rasters)
- **lattice**: Class-I triangular; F^2 triangles/face; lower [i,j]→[i+1,j]→[i,j+1], upper [i+1,j]→[i+1,j+1]→[i,j+1]
- **bary→lattice rule**: u=αF, v=βF; upper if frac(u)+frac(v) > 1−ε; clamp if floor(u)+floor(v) ≥ F

Face 0 neighbors and permutations:
- **BC (opp A)**: neighbor=1, perm=(A,B,C)->(B',A',C')
- **CA (opp B)**: neighbor=4, perm=(A,B,C)->(B',A',C')
- **AB (opp C)**: neighbor=6, perm=(A,B,C)->(A',C',B')

Short table for F=4, face 0:

```
corners: [[-0.850651,0,0.525731],[0,0.525731,0.850651],[-0.525731,0.850651,0]]
first_10_face_vertices: [{"i":0,"j":0,"id":11},{"i":0,"j":1,"id":16},{"i":0,"j":2,"id":20},{"i":0,"j":3,"id":23},{"i":0,"j":4,"id":5},{"i":1,"j":0,"id":14},{"i":1,"j":1,"id":18},{"i":1,"j":2,"id":21},{"i":1,"j":3,"id":22},{"i":2,"j":0,"id":13}]
first_8_triangles: [{"kind":"lower","idx":[16,14,11]},{"kind":"upper","idx":[16,18,14]},{"kind":"lower","idx":[20,18,16]},{"kind":"upper","idx":[20,21,18]},{"kind":"lower","idx":[23,21,20]},{"kind":"upper","idx":[23,22,21]},{"kind":"lower","idx":[5,22,23]},{"kind":"lower","idx":[18,13,14]}]
```

Assumptions/notes:
- **F power-of-two**: THREE.IcosahedronGeometry(detail) implies F=2^detail.
- **Diagonal tie**: uses frac(u)+frac(v) > 1−ε for upper; else lower.
- **Neighbor perms**: derived by shared-vertex correspondence per base icosa faces.