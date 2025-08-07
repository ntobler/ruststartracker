pub struct OrderedCombinations<const K: usize> {
    n: u32,
    indices: Option<[u32; K]>,
}

impl<const K: usize> OrderedCombinations<K> {
    pub fn new(n: u32) -> Result<Self, &'static str> {
        if n < K as u32 {
            return Err("n must be at least K");
        }

        let mut initial = [0u32; K];
        for i in 0..K {
            initial[i] = i as u32;
        }

        Ok(Self {
            n,
            indices: Some(initial),
        })
    }
}

impl<const K: usize> Iterator for OrderedCombinations<K> {
    type Item = [u32; K];

    fn next(&mut self) -> Option<Self::Item> {
        let ret = self.indices?;
        let mut indices = self.indices?;

        // Increment indices
        for k in (0..K - 1).rev() {
            if indices[k] != indices[K - 1] - (K - 1 - k) as u32 {
                indices[k] += 1;
                for k2 in k + 1..K - 1 {
                    indices[k2] = indices[k2 - 1] + 1;
                }
                self.indices = Some(indices);
                return Some(ret);
            }
        }

        // last position is special as it affect the max index
        if indices[K - 1] != self.n - 1 {
            // pos 012 limit, increase max, reinit [0..k-2, max]
            for i in 0..K - 1 {
                indices[i] = i as u32;
            }
            indices[K - 1] += 1;

            self.indices = Some(indices);
            return Some(ret);
        }

        self.indices = None;
        return Some(ret);
    }
}

#[cfg(test)]
mod tests {

    use itertools::Itertools;

    use super::*;

    #[test]
    fn test_combination_order_old() {
        let n = 5;
        let iter = itertools::Itertools::combinations(0..n, 3).map(|item| item.try_into().unwrap());
        let v: Vec<[u32; 3]> = iter.collect_vec();
        assert_eq!(v[0], [0, 1, 2]);
        assert_eq!(v[1], [0, 1, 3]);
        assert_eq!(v[2], [0, 1, 4]);
        assert_eq!(v[3], [0, 2, 3]);
        assert_eq!(v[4], [0, 2, 4]);
        assert_eq!(v[5], [0, 3, 4]);
        assert_eq!(v[6], [1, 2, 3]);
        assert_eq!(v[7], [1, 2, 4]);
        assert_eq!(v[8], [1, 3, 4]);
        assert_eq!(v[9], [2, 3, 4]);
        assert_eq!(v.len(), 10);
    }

    #[test]
    fn test_combination_order_sfa_3_5() {
        let n = 5;
        let iter = OrderedCombinations::<3>::new(n).unwrap();
        let v: Vec<[u32; 3]> = iter.take(11).collect();

        assert_eq!(v[0], [0, 1, 2]);
        assert_eq!(v[1], [0, 1, 3]);
        assert_eq!(v[2], [0, 2, 3]);
        assert_eq!(v[3], [1, 2, 3]);
        assert_eq!(v[4], [0, 1, 4]);
        assert_eq!(v[5], [0, 2, 4]);
        assert_eq!(v[6], [0, 3, 4]);
        assert_eq!(v[7], [1, 2, 4]);
        assert_eq!(v[8], [1, 3, 4]);
        assert_eq!(v[9], [2, 3, 4]);
        assert_eq!(v.len(), 10);
    }

    #[test]
    fn test_combination_order_3_6() {
        let n = 6;
        let iter = OrderedCombinations::<3>::new(n).unwrap();
        let v: Vec<[u32; 3]> = iter.take(21).collect();

        assert_eq!(v[0], [0, 1, 2]);
        assert_eq!(v[1], [0, 1, 3]);
        assert_eq!(v[2], [0, 2, 3]);
        assert_eq!(v[3], [1, 2, 3]);
        assert_eq!(v[4], [0, 1, 4]);
        assert_eq!(v[5], [0, 2, 4]);
        assert_eq!(v[6], [0, 3, 4]);
        assert_eq!(v[7], [1, 2, 4]);
        assert_eq!(v[8], [1, 3, 4]);
        assert_eq!(v[9], [2, 3, 4]);
        assert_eq!(v[10], [0, 1, 5]);
        assert_eq!(v[11], [0, 2, 5]);
        assert_eq!(v[12], [0, 3, 5]);
        assert_eq!(v[13], [0, 4, 5]);
        assert_eq!(v[14], [1, 2, 5]);
        assert_eq!(v[15], [1, 3, 5]);
        assert_eq!(v[16], [1, 4, 5]);
        assert_eq!(v[17], [2, 3, 5]);
        assert_eq!(v[18], [2, 4, 5]);
        assert_eq!(v[19], [3, 4, 5]);
        assert_eq!(v.len(), 20);
    }

    #[test]
    fn test_combination_order_4_5() {
        let n = 5;
        let iter = OrderedCombinations::<4>::new(n).unwrap();
        let v: Vec<[u32; 4]> = iter.take(12).collect();

        assert_eq!(v[0], [0, 1, 2, 3]);
        assert_eq!(v[1], [0, 1, 2, 4]);
        assert_eq!(v[2], [0, 1, 3, 4]);
        assert_eq!(v[3], [0, 2, 3, 4]);
        assert_eq!(v[4], [1, 2, 3, 4]);
        assert_eq!(v.len(), 5);
    }

    #[test]
    fn test_combination_order_4_6() {
        let n = 6;
        let iter = OrderedCombinations::<4>::new(n).unwrap();
        let v: Vec<[u32; 4]> = iter.collect();

        assert_eq!(v[0], [0, 1, 2, 3]);
        assert_eq!(v[1], [0, 1, 2, 4]);
        assert_eq!(v[2], [0, 1, 3, 4]);
        assert_eq!(v[3], [0, 2, 3, 4]);
        assert_eq!(v[4], [1, 2, 3, 4]);
        assert_eq!(v[5], [0, 1, 2, 5]);
        assert_eq!(v[6], [0, 1, 3, 5]);
        assert_eq!(v[7], [0, 1, 4, 5]);
        assert_eq!(v[8], [0, 2, 3, 5]);
        assert_eq!(v[9], [0, 2, 4, 5]);
        assert_eq!(v[10], [0, 3, 4, 5]);
        assert_eq!(v[11], [1, 2, 3, 5]);
        assert_eq!(v[12], [1, 2, 4, 5]);
        assert_eq!(v[13], [1, 3, 4, 5]);
        assert_eq!(v[14], [2, 3, 4, 5]);
        assert_eq!(v.len(), 15);
    }

    #[test]
    fn test_combination_order_3_3() {
        let iter = OrderedCombinations::<3>::new(3).unwrap();
        let v: Vec<[u32; 3]> = iter.collect();
        assert_eq!(v[0], [0, 1, 2]);
        assert_eq!(v.len(), 1);
    }

    #[test]
    fn test_combination_order_3_2() {
        assert!(OrderedCombinations::<3>::new(2).is_err());
    }
}
