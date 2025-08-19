#[cfg(test)]
mod tests {
    use super::*;
    use crate::icosa::to_gpu_faces;
    use crate::math::*;

    #[test]
    fn face_table_is_deterministic() {
        let a = crate::build_face_table();
        let b = crate::build_face_table();
        assert_eq!(a.len(), 20);
        for i in 0..20 {
            assert!((a[i].n.x - b[i].n.x).abs() < 1e-8);
        }
        let (_fa, na) = to_gpu_faces(&a);
        let (_fb, nb) = to_gpu_faces(&b);
        assert_eq!(na.len(), nb.len());
    }

    #[test]
    fn normals_point_outward() {
        let faces = crate::build_face_table();
        for f in faces {
            let centroid = f.a.add(f.b).add(f.c);
            assert!(f.n.dot(centroid) > 0.0);
        }
    }

    #[test]
    fn neighbor_slots_are_filled() {
        let faces = crate::build_face_table();
        for f in faces {
            for k in 0..3 {
                assert_ne!(f.neighbor_opp[k], u32::MAX);
            }
        }
    }
}


