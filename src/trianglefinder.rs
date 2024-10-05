pub struct TriangleFinder {
    connections_ab: Vec<[u32; 2]>,
    connections_ac: Vec<[u32; 2]>,
    connections_bc: Vec<[u32; 2]>,
}

impl TriangleFinder {
    pub fn new(
        connections_ab: Vec<[u32; 2]>,
        connections_ac: Vec<[u32; 2]>,
        connections_bc: Vec<[u32; 2]>,
    ) -> Self {
        TriangleFinder {
            connections_ab,
            connections_ac,
            connections_bc,
        }
    }

    pub fn get(&self) -> Option<[u32; 3]> {
        let mut c_candidates = Vec::new();
        let mut flipped_c_candidates = Vec::new();

        for [a, b] in &self.connections_ab {
            c_candidates.clear();
            flipped_c_candidates.clear();

            for [n1, n2] in &self.connections_ac {
                if a == n1 {
                    c_candidates.push(n2)
                } else if a == n2 {
                    c_candidates.push(n1)
                }

                if b == n1 {
                    flipped_c_candidates.push(n2)
                } else if b == n2 {
                    flipped_c_candidates.push(n1)
                }
            }

            if c_candidates.is_empty() && flipped_c_candidates.is_empty() {
                continue;
            }

            for [n1, n2] in &self.connections_bc {
                if b == n1 {
                    c_candidates.push(n2)
                } else if b == n2 {
                    c_candidates.push(n1)
                }

                if a == n1 {
                    flipped_c_candidates.push(n2)
                } else if a == n2 {
                    flipped_c_candidates.push(n1)
                }
            }

            if c_candidates.len() >= 2 {
                c_candidates.sort_unstable();
                let mut prev = c_candidates[0];
                for c in c_candidates.iter().skip(1) {
                    if *c == prev {
                        return Some([*a, *b, **c]);
                    }
                    prev = *c;
                }
            }

            if flipped_c_candidates.len() >= 2 {
                flipped_c_candidates.sort_unstable();
                let mut prev = flipped_c_candidates[0];
                for c in flipped_c_candidates.iter().skip(1) {
                    if *c == prev {
                        return Some([*a, *b, **c]);
                    }
                    prev = *c;
                }
            }
        }
        return None;
    }

    pub fn get_all(&self) -> Vec<[u32; 3]> {
        let mut c_candidates = Vec::new();
        let mut flipped_c_candidates = Vec::new();
        let mut results = Vec::new();

        for [a, b] in &self.connections_ab {
            c_candidates.clear();
            flipped_c_candidates.clear();

            for [n1, n2] in &self.connections_ac {
                if a == n1 {
                    c_candidates.push(n2)
                } else if a == n2 {
                    c_candidates.push(n1)
                }

                if b == n1 {
                    flipped_c_candidates.push(n2)
                } else if b == n2 {
                    flipped_c_candidates.push(n1)
                }
            }

            if c_candidates.is_empty() && flipped_c_candidates.is_empty() {
                continue;
            }

            for [n1, n2] in &self.connections_bc {
                if b == n1 {
                    c_candidates.push(n2)
                } else if b == n2 {
                    c_candidates.push(n1)
                }

                if a == n1 {
                    flipped_c_candidates.push(n2)
                } else if a == n2 {
                    flipped_c_candidates.push(n1)
                }
            }

            if c_candidates.len() >= 2 {
                c_candidates.sort_unstable();
                let mut prev = c_candidates[0];
                for c in c_candidates.iter().skip(1) {
                    if *c == prev {
                        results.push([*a, *b, **c]);
                    }
                    prev = *c;
                }
            }

            if flipped_c_candidates.len() >= 2 {
                flipped_c_candidates.sort_unstable();
                let mut prev = flipped_c_candidates[0];
                for c in flipped_c_candidates.iter().skip(1) {
                    if *c == prev {
                        results.push([*a, *b, **c]);
                    }
                    prev = *c;
                }
            }
        }
        results
    }
}

pub struct IterTriangleFinder {
    triangle_finder: TriangleFinder,
    segment_id: usize,
    c_candidates: Vec<u32>,
    flipped_c_candidates: Vec<u32>,
    result: Vec<[u32; 3]>,
    index: usize,
}

impl IterTriangleFinder {
    pub fn new(triangle_finder: TriangleFinder) -> Self {
        IterTriangleFinder {
            triangle_finder,
            segment_id: 0,
            c_candidates: Vec::new(),
            flipped_c_candidates: Vec::new(),
            result: Vec::new(),
            index: 0,
        }
    }

    fn search_next(&mut self) -> bool {
        let finder = &self.triangle_finder;

        self.result.clear();

        while self.segment_id < finder.connections_ab.len() {
            let [a, b] = finder.connections_ab[self.segment_id];

            self.segment_id += 1;

            self.c_candidates.clear();
            self.flipped_c_candidates.clear();

            for [n1, n2] in finder.connections_ac.iter() {
                if a == *n1 {
                    self.c_candidates.push(*n2)
                } else if a == *n2 {
                    self.c_candidates.push(*n1)
                }

                if b == *n1 {
                    self.flipped_c_candidates.push(*n2)
                } else if b == *n2 {
                    self.flipped_c_candidates.push(*n1)
                }
            }

            if self.c_candidates.is_empty() && self.flipped_c_candidates.is_empty() {
                continue;
            }

            if self.segment_id >= finder.connections_ab.len() {
                return false;
            }

            for [n1, n2] in finder.connections_bc.iter() {
                if b == *n1 {
                    self.c_candidates.push(*n2)
                } else if b == *n2 {
                    self.c_candidates.push(*n1)
                }

                if a == *n1 {
                    self.flipped_c_candidates.push(*n2)
                } else if a == *n2 {
                    self.flipped_c_candidates.push(*n1)
                }
            }

            if self.c_candidates.len() >= 2 {
                self.c_candidates.sort_unstable();
                let mut old = self.c_candidates[0];
                for c in self.c_candidates.iter().skip(1) {
                    if *c == old {
                        self.result.push([a, b, *c]);
                    }
                    old = *c;
                }
            }

            if self.flipped_c_candidates.len() >= 2 {
                self.flipped_c_candidates.sort_unstable();
                let mut old = self.flipped_c_candidates[0];
                for c in self.flipped_c_candidates.iter().skip(1) {
                    if *c == old {
                        self.result.push([a, b, *c]);
                    }
                    old = *c;
                }
            }

            if !self.result.is_empty() {
                return true;
            }
        }
        return false;
    }
}

impl Iterator for IterTriangleFinder {
    type Item = [u32; 3];
    fn next(&mut self) -> Option<[u32; 3]> {
        loop {
            if self.index >= self.result.len() {
                self.index = 0;
                if !self.search_next() {
                    return None;
                }
            } else {
                let res = Some(self.result[self.index]);
                self.index += 1;
                return res;
            }
        }
    }
}
