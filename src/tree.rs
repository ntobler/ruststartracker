use core::f32;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;

pub struct UnitVectorLookup {
    kdtree: KdTree<f32, usize, [f32; 3]>,
}

pub fn dot_product(a: &[f32], b: &[f32]) -> f32 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

impl UnitVectorLookup {
    pub fn new(vectors: &[[f32; 3]]) -> UnitVectorLookup {
        let mut kdtree = KdTree::new(3);
        for (i, item) in vectors.iter().enumerate() {
            kdtree.add(*item, i).unwrap();
        }
        UnitVectorLookup { kdtree }
    }

    pub fn lookup_nearest(&self, vector: &[f32; 3]) -> usize {
        let res = self.kdtree.nearest(vector, 1, &squared_euclidean).unwrap();
        *(res[0].1)
    }

    pub fn look_up_close_angles(
        &self,
        vectors: &[[f32; 3]],
        max_angle_rad: f32,
    ) -> Vec<([u32; 2], f32)> {
        let threshold = maths_rs::cos(max_angle_rad);
        let mut index_pairs = Vec::new();
        for a in 0..vectors.len() {
            let vec_a = &vectors[a];

            for (_, b) in self.kdtree.iter_nearest(vec_a, &squared_euclidean).unwrap() {
                let vec_b = &vectors[*b];
                let dotp = dot_product(vec_a, vec_b);
                if dotp < threshold {
                    break;
                }
                if a < *b {
                    index_pairs
                        .push(([a as u32, *b as u32], maths_rs::acos(dotp.clamp(-1.0, 1.0))));
                }
            }
        }
        index_pairs
    }
}
