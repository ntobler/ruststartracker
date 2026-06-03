use core::f32;

use kdtree::distance::squared_euclidean;
use kdtree::KdTree;

pub struct UnitVectorLookup {
    kdtree: KdTree<f32, usize, [f32; 3]>,
}

#[inline(always)]
pub fn dot_product(a: &[f32; 3], b: &[f32; 3]) -> f32 {
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
        // kdtree.nearest is equally fast as kdtree.iter_nearest
        let res = self.kdtree.nearest(vector, 1, &squared_euclidean).unwrap();
        *(res[0].1)
    }

    pub fn look_up_close_angles(
        &self,
        vectors: &[[f32; 3]],
        magnitudes: &[f32],
        cos_max_angle: f32,
        max_magnitude: f32,
    ) -> Vec<([u32; 2], f32)> {
        let mut index_pairs = Vec::new();
        for a in 0..vectors.len() {
            if magnitudes[a] > max_magnitude {
                continue;
            }
            let vec_a = &vectors[a];
            for (_, b) in self.kdtree.iter_nearest(vec_a, &squared_euclidean).unwrap() {
                if magnitudes[*b] > max_magnitude {
                    continue;
                }
                let vec_b = &vectors[*b];
                let dotp = dot_product(vec_a, vec_b);
                if dotp < cos_max_angle {
                    // If angle is too large, break here
                    break;
                }
                if a < *b {
                    index_pairs.push(([a as u32, *b as u32], dotp.clamp(-1.0, 1.0)));
                }
            }
        }
        index_pairs
    }
}

#[cfg(test)]
mod tests {

    use rand::rng;
    use rand_distr::{Distribution, Normal};

    use super::*;

    #[test]
    fn test_tree() {
        let mut rng = rng();
        let normal = Normal::new(0.0, 1.0).unwrap(); // mean = 0, std dev = 1

        let samples: Vec<[f32; 3]> = (0..100)
            .map(|_| {
                let x = normal.sample(&mut rng) as f32;
                let y = normal.sample(&mut rng) as f32;
                let z = normal.sample(&mut rng) as f32;
                let mag = f32::sqrt(x * x + y * y + z * z);
                [x / mag, y / mag, z / mag]
            })
            .collect();

        let lookup = UnitVectorLookup::new(&samples);

        assert!(lookup.lookup_nearest(&samples[0]) == 0);
    }
}
