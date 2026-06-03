pub struct IterTriangleFinder<'a> {
    connections_ab: &'a [[u32; 2]],
    connections_ac: &'a [[u32; 2]],
    connections_bc: &'a [[u32; 2]],
    ab_index: usize,
    c_candidates_given_a_is_ab1: Vec<u32>,
    c_candidates_given_a_is_ab2: Vec<u32>,
    matches: Vec<[u32; 3]>,
}

impl<'a> IterTriangleFinder<'a> {
    pub fn new(
        connections_ab: &'a [[u32; 2]],
        connections_ac: &'a [[u32; 2]],
        connections_bc: &'a [[u32; 2]],
    ) -> Self {
        IterTriangleFinder {
            connections_ab,
            connections_ac,
            connections_bc,
            ab_index: 0,
            c_candidates_given_a_is_ab1: Vec::with_capacity(128),
            c_candidates_given_a_is_ab2: Vec::with_capacity(128),
            matches: Vec::with_capacity(128),
        }
    }

    fn search_next(&mut self, ab_index: usize) {
        // We are looking for triples (a, b, c) such that:
        // (a, b) in connections_ab
        // (a, c) in connections_ac
        // (b, c) in connections_bc
        //
        // Creating an index for connections_ac and connections_bc would be faster,
        // but it would also require more memory. Instead, we will just iterate over
        // the connections_ac and connections_bc for each (a, b) pair in connections_ab.

        let [ab1, ab2] = self.connections_ab[ab_index];

        self.c_candidates_given_a_is_ab1.clear();
        self.c_candidates_given_a_is_ab2.clear();

        for &[ac1, ac2] in self.connections_ac {
            if ab1 == ac1 {
                self.c_candidates_given_a_is_ab1.push(ac2)
            } else if ab1 == ac2 {
                self.c_candidates_given_a_is_ab1.push(ac1)
            }

            if ab2 == ac1 {
                self.c_candidates_given_a_is_ab2.push(ac2)
            } else if ab2 == ac2 {
                self.c_candidates_given_a_is_ab2.push(ac1)
            }
        }

        for &[bc1, bc2] in self.connections_bc {
            if ab2 == bc1 {
                if self.c_candidates_given_a_is_ab1.contains(&bc2) {
                    self.matches.push([ab1, ab2, bc2]);
                }
            } else if ab2 == bc2 {
                if self.c_candidates_given_a_is_ab1.contains(&bc1) {
                    self.matches.push([ab1, ab2, bc1]);
                }
            } else if ab1 == bc1 {
                if self.c_candidates_given_a_is_ab2.contains(&bc2) {
                    self.matches.push([ab2, ab1, bc2]);
                }
            } else if ab1 == bc2 {
                if self.c_candidates_given_a_is_ab2.contains(&bc1) {
                    self.matches.push([ab2, ab1, bc1]);
                }
            }
        }
    }
}

impl<'a> Iterator for IterTriangleFinder<'a> {
    type Item = [u32; 3];
    fn next(&mut self) -> Option<[u32; 3]> {
        loop {
            match self.matches.pop() {
                Some(res) => return Some(res),
                None => {
                    if self.ab_index < self.connections_ab.len() {
                        self.search_next(self.ab_index);
                        self.ab_index += 1;
                    } else {
                        return None;
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_triangle_finder() {
        let ab = vec![[234, 5643], [1, 2], [2, 4], [3, 9], [2, 6], [1, 0], [3, 1]];
        let ac = vec![[345, 2343], [8, 2], [3, 4], [1, 7], [0, 5], [3, 1]];
        let bc = vec![[435, 4355], [1, 0], [4, 8], [8, 1], [1, 9], [2, 6], [3, 7]];

        let i = IterTriangleFinder::new(&ab, &ac, &bc);
        let vec = i.collect::<Vec<[u32; 3]>>();
        assert!(vec == vec![[2, 1, 8], [2, 4, 8], [3, 9, 1], [1, 3, 7]]);
    }

    #[test]
    fn test_triangle_finder_performance() {
        use rand::{rngs::StdRng, Rng, SeedableRng};

        const N: usize = 200;
        const M: u32 = 200;

        let mut rng = StdRng::seed_from_u64(42);
        let ab: Vec<[u32; 2]> = (0..N)
            .map(|_| [rng.random_range(0..M), rng.random_range(0..M)])
            .filter(|&[a, b]| a != b)
            .collect();
        let ac: Vec<[u32; 2]> = (0..N)
            .map(|_| [rng.random_range(0..M), rng.random_range(0..M)])
            .filter(|&[a, b]| a != b)
            .collect();
        let bc: Vec<[u32; 2]> = (0..N)
            .map(|_| [rng.random_range(0..M), rng.random_range(0..M)])
            .filter(|&[a, b]| a != b)
            .collect();

        let i = IterTriangleFinder::new(&ab, &ac, &bc);

        let start = std::time::Instant::now();
        let vec = i.collect::<Vec<[u32; 3]>>();
        let duration = start.elapsed();
        println!("Vec length: {}", vec.len());
        println!("Collection took: {:?}", duration);
    }
}
