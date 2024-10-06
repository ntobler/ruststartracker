use crate::tree;
use crate::trianglefinder;
use std::time::Instant;

extern crate nalgebra as na;

pub struct StarMatcher {
    stars_xyz: Vec<[f32; 3]>,
    star_index: tree::UnitVectorLookup,
    inter_star_angles: Vec<f32>,
    inter_star_angle_pairs: Vec<[u32; 2]>,
    /// tolerance of inter star angle in rad
    inter_star_angle_tolerance: f32,
    // polynom with terms [c0, c1, c2, ..] (c0 + x c1 + x^2 c2)
    inter_star_index_polygon: Vec<f32>,
    max_inter_star_angle: f32,
    n_minimum_matches: usize,
    timeout_secs: f32,
}

/// Return angle between two normalized 3-dimensional vectors.
fn angle(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    maths_rs::acos(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
}

/// Evaluate a polynom
fn polyval(coeffs: &[f32], x: f32) -> f32 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

/// Naive implementation for tree::UnitVectorLookup::look_up_close_angles
/// Calculates all star angles, instead of looking up neighbors in a spatial index
/// Turned out to be faster in some cases
pub fn look_up_close_angles_naive(
    vectors: &[[f32; 3]],
    max_angle_rad: f32,
) -> Vec<([u32; 2], f32)> {
    let threshold = maths_rs::cos(max_angle_rad);
    let mut index_pairs = Vec::new();
    for a in 0..vectors.len() {
        let vec_a = &vectors[a];
        for b in (a + 1)..vectors.len() {
            let vec_b = &vectors[b];
            let dotp = tree::dot_product(vec_a, vec_b);
            if dotp >= threshold {
                index_pairs.push(([a as u32, b as u32], maths_rs::acos(dotp)));
            }
        }
    }
    index_pairs
}

pub fn get_inter_star_index(
    star_index: &tree::UnitVectorLookup,
    stars_xyz: &[[f32; 3]],
    max_angle_rad: f32,
) -> (Vec<[u32; 2]>, Vec<f32>, Vec<f32>) {
    let mut index_pairs = star_index.look_up_close_angles(stars_xyz, max_angle_rad);

    index_pairs.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    let angles: Vec<f32> = index_pairs.iter().map(|x| x.1).collect();
    let indices: Vec<f32> = (0..angles.len()).map(|i| i as f32).collect();

    // polynom with terms [c0, c1, c2, ..] (c0 + x c1 + x^2 c2)
    let mut polynom = polyfit_rs::polyfit_rs::polyfit(&angles, &indices, 2).unwrap();

    let errors: Vec<f32> = angles
        .iter()
        .enumerate()
        .map(|x| polyval(&polynom, *x.1) - (x.0 as f32))
        .collect();

    // let min = errors.iter().cloned().fold(f32::INFINITY, f32::min);
    let max = errors.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

    polynom[0] -= max;

    let pairs = index_pairs.iter().map(|x| x.0).collect();

    (pairs, angles, polynom)
}

impl StarMatcher {
    pub fn new(
        stars_xyz: Vec<[f32; 3]>,
        max_inter_star_angle: f32,
        inter_star_angle_tolerance: f32,
        n_minimum_matches: usize,
        timeout_secs: f32,
    ) -> Self {
        let star_index = tree::UnitVectorLookup::new(&stars_xyz);

        let (pairs, angles, polynom) =
            get_inter_star_index(&star_index, &stars_xyz, max_inter_star_angle);

        StarMatcher {
            stars_xyz,
            star_index,
            inter_star_angles: angles,
            inter_star_angle_pairs: pairs,
            inter_star_angle_tolerance,
            inter_star_index_polygon: polynom,
            max_inter_star_angle,
            n_minimum_matches,
            timeout_secs,
        }
    }

    pub fn find(&self, obs_xyz: Vec<[f32; 3]>) -> Result<MatchResult, &'static str> {
        let start_instant = Instant::now();

        for obs_indices in self.triangle_combinations_iterator(obs_xyz.len() as u32) {
            let [a, b, c] = obs_indices;

            // Get positions of the 3 observations
            let obs_a = obs_xyz[a as usize];
            let obs_b = obs_xyz[b as usize];
            let obs_c = obs_xyz[c as usize];

            // find the inter star angles between the 3 observations
            let angle_ab = angle(&obs_a, &obs_b);
            let angle_ac = angle(&obs_a, &obs_c);
            let angle_bc = angle(&obs_b, &obs_c);

            // Continue if any angle is larger than the maximum
            // angle contained in the inter star angle index
            let t = self.max_inter_star_angle;
            if (angle_ab > t) || (angle_ac > t) || (angle_bc > t) {
                continue;
            }

            // Look up all pairs in the catalog that could match the observations
            let ab_pairs = self.pair_lookup(angle_ab);
            let ac_pairs = self.pair_lookup(angle_ac);
            let bc_pairs = self.pair_lookup(angle_bc);

            let finder = trianglefinder::TriangleFinder::new(
                ab_pairs.to_vec(),
                ac_pairs.to_vec(),
                bc_pairs.to_vec(),
            );

            let iter_finder = trianglefinder::IterTriangleFinder::new(finder);

            // Iterate over possible matching triangles
            for value in iter_finder {
                match self.check(&obs_indices, &obs_xyz, &value) {
                    None => {}
                    Some(x) => return Ok(x),
                };
            }

            if start_instant.elapsed().as_secs_f32() > self.timeout_secs {
                return Err("Timeout reached");
            }
        }
        Err("Search exhausted")
    }

    /// Iterator over combinations of stars forming triangles.
    fn triangle_combinations_iterator(&self, n: u32) -> impl Iterator<Item = [u32; 3]> {
        itertools::Itertools::combinations(0..n, 3).map(|item| item.try_into().unwrap())
    }

    /// Get star pairs that match given inter star angle.
    fn pair_lookup(&self, inter_star_angle: f32) -> &[[u32; 2]] {
        let lower_threshold =
            maths_rs::max(inter_star_angle - self.inter_star_angle_tolerance, 0.0);
        let upper_threshold =
            maths_rs::max(inter_star_angle + self.inter_star_angle_tolerance, 0.0);
        let lower_index_float = polyval(&self.inter_star_index_polygon, lower_threshold);
        let upper_index_float = polyval(&self.inter_star_index_polygon, upper_threshold);
        let max = self.inter_star_angles.len() - 1;
        let mut lower_index = (lower_index_float as usize).clamp(0, max);
        let mut upper_index = (upper_index_float as usize).clamp(0, max);

        lower_index = maths_rs::min(lower_index, upper_index);

        // print!("Pair look up of angle {:?}\n", inter_star_angle);
        // let lower_prev = lower_index;
        while lower_index < max {
            if self.inter_star_angles[lower_index] > lower_threshold {
                break;
            } else {
                lower_index += 1;
            }
        }
        // print!(
        //     "found lower index {:?} in {:?} steps, got {:?}rad when looking for {:?}rad\n",
        //     lower_index,
        //     lower_index - lower_prev,
        //     self.inter_star_angles[lower_index],
        //     lower_threshold
        // );
        // let upper_prev = upper_index;
        while upper_index < max {
            if self.inter_star_angles[upper_index] > upper_threshold {
                break;
            } else {
                upper_index += 1;
            }
        }
        // print!(
        //     "found lower index {:?} in {:?} steps, got {:?}rad when looking for {:?}rad\n",
        //     upper_index,
        //     upper_index - upper_prev,
        //     self.inter_star_angles[upper_index],
        //     upper_threshold
        // );

        &self.inter_star_angle_pairs[lower_index..upper_index]
    }

    fn check(
        &self,
        obs_indices: &[u32; 3],
        obs_xyz: &[[f32; 3]],
        cat_indices: &[u32; 3],
    ) -> Option<MatchResult> {
        // Get vectors of observed triangle
        let obs_triangle_xyz = [
            obs_xyz[obs_indices[0] as usize],
            obs_xyz[obs_indices[1] as usize],
            obs_xyz[obs_indices[2] as usize],
        ];

        // Get vectors of matched catalog triangle
        let cat_triangle_xyz = [
            self.stars_xyz[cat_indices[0] as usize],
            self.stars_xyz[cat_indices[1] as usize],
            self.stars_xyz[cat_indices[2] as usize],
        ];

        // Fit rotation matrix on triangle
        let rotm = match attitude_svd(&cat_triangle_xyz, &obs_triangle_xyz) {
            None => {
                print!("1st Attitude svd failed\n");
                return None;
            }
            Some(value) => value.cast::<f32>(),
        };

        // Transform all observations, such that they align with the catalog stars
        let obs_xyz_mat =
            na::Matrix3xX::from_iterator(obs_xyz.len(), obs_xyz.into_iter().flatten().cloned());
        let obs_transformed = rotm * obs_xyz_mat;

        // Find close neighbors for each observation
        let dotp_threshold = maths_rs::cos(self.inter_star_angle_tolerance);
        let mut selected_obs_xyz = Vec::new();
        let mut selected_cat_xyz = Vec::new();
        let mut selected_cat_indices = Vec::new();
        for obs_i in 0..obs_transformed.ncols() {
            // version 1: use spatial look up
            let obs_vec = obs_transformed.column(obs_i);
            let obs = [obs_vec[0], obs_vec[1], obs_vec[2]];

            //Look up closest star in the catalog to the transformed position of the observation
            let closest_index = self.star_index.lookup_nearest(&obs);

            // Use star if it is close than the allowed threshold
            let closest_cat_star = &self.stars_xyz[closest_index];
            if tree::dot_product(&obs, closest_cat_star) >= dotp_threshold {
                selected_cat_xyz.push(*closest_cat_star);
                selected_obs_xyz.push(obs_xyz[obs_i]);
                selected_cat_indices.push(closest_index as u32)
            }
        }

        // Do not proceed if there are less than the minimum required amount of stars
        if selected_cat_xyz.len() < self.n_minimum_matches {
            // print!("Less than {} close neighbors found\n", self.n_minimum_matches);
            return None;
        }

        // Fit rotation matrix on selected observations
        let final_rotm = match attitude_svd(&selected_cat_xyz, &selected_obs_xyz) {
            None => {
                print!("2nd Attitude svd failed\n");
                return None;
            }
            Some(value) => value,
        };

        // Convert rotation matrix into quaternion
        let quat = na::UnitQuaternion::from_rotation_matrix(&na::Rotation3::from_matrix_unchecked(
            final_rotm,
        ));

        //Return result with statistics
        Some(MatchResult {
            quat: [quat.i as f32, quat.j as f32, quat.k as f32, quat.w as f32],
            match_ids: selected_cat_indices,
            n_matches: selected_cat_xyz.len() as u32,
            obs_matched: selected_cat_xyz,
        })
    }
}

pub struct MatchResult {
    pub quat: [f32; 4],
    pub match_ids: Vec<u32>,
    pub n_matches: u32,
    pub obs_matched: Vec<[f32; 3]>,
}

/// Solve Wahba's problem using SVD method.
///
/// Rotation matrix from inertial frame to sensor frame
///
/// G. Wahba, "A Least Squares Estimate of Spacecraft Attitude", (1965)
///
/// # Arguments
///
/// - `cat_xyz` - Catalog unit vectors in the inertial frame, 3xn matrix
/// - `obs_xyz``- Measured unit vectors in the sensor frame, 3xn matrix
pub fn attitude_svd(cat_xyz: &[[f32; 3]], obs_xyz: &[[f32; 3]]) -> Option<na::Matrix3<f64>> {
    let mut mat = na::Matrix3::<f64>::zeros();

    for i in 0..cat_xyz.len() {
        let [d, e, f] = cat_xyz[i];
        let [a, b, c] = obs_xyz[i];
        let outer_prod = na::Matrix3::new(
            a * d,
            b * d,
            c * d,
            a * e,
            b * e,
            c * e,
            a * f,
            b * f,
            c * f,
        );
        mat += outer_prod.try_cast::<f64>().unwrap();
    }
    // Perform SVD
    let svd = mat.svd(true, true);

    let u = svd.u.unwrap();
    let v = svd.v_t.unwrap();
    let d = u.determinant() * v.determinant();
    let m = na::Matrix3::from_diagonal(&na::Vector3::new(1.0, 1.0, d));
    let t = (u * m) * v;
    match t.determinant() >= 0.0 {
        true => Some(t),
        false => None,
    }
}
