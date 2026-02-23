use std::fmt;

use crate::tree;
use crate::trianglefinder;
use std::time::Instant;

extern crate nalgebra as na;

/// Return angle between two normalized 3-dimensional vectors.
fn angle(a: &[f32; 3], b: &[f32; 3]) -> f32 {
    maths_rs::acos(a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
}

/// Evaluate a polynomial
fn polyval(coeffs: &[f32], x: f32) -> f32 {
    coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
}

/// Naive implementation for tree::UnitVectorLookup::look_up_close_angles
/// Calculates all star angles, instead of looking up neighbors in a spatial index
/// Turned out to be faster in some cases
pub fn look_up_close_angles_naive(
    vectors: &[[f32; 3]],
    magnitudes: &[f32],
    max_angle_rad: f32,
    max_magnitude: f32,
) -> Vec<([u32; 2], f32)> {
    let threshold = maths_rs::cos(max_angle_rad);
    let mut index_pairs = Vec::new();
    for a in 0..vectors.len() {
        if magnitudes[a] > max_magnitude {
            continue;
        }
        let vec_a = &vectors[a];
        for b in (a + 1)..vectors.len() {
            if magnitudes[b] > max_magnitude {
                continue;
            }
            let vec_b = &vectors[b];
            let dotp = tree::dot_product(vec_a, vec_b);
            if dotp >= threshold {
                index_pairs.push(([a as u32, b as u32], maths_rs::acos(dotp)));
            }
        }
    }
    index_pairs
}

pub struct InterStarIndex {
    pub pairs: Vec<[u32; 2]>,
    pub angles: Vec<f32>,
    // polynomial with terms [c0, c1, c2, ..] (c0 + x c1 + x^2 c2)
    pub polynomial: [f32; 3],
}

impl InterStarIndex {
    pub fn new(
        star_index: &tree::UnitVectorLookup,
        stars_xyz: &[[f32; 3]],
        stars_mag: &[f32],
        max_angle_rad: f32,
        max_magnitude: f32,
    ) -> Result<InterStarIndex, &'static str> {
        if stars_xyz.len() != stars_mag.len() {
            return Err("stars_xyz and stars_mag must have same length");
        }

        let mut index_pairs =
            star_index.look_up_close_angles(stars_xyz, stars_mag, max_angle_rad, max_magnitude);

        if index_pairs.is_empty() {
            return Err("Given star positions do not result in any angles below the threshold.");
        }

        index_pairs.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        let angles: Vec<f32> = index_pairs.iter().map(|x| x.1).collect();
        let indices: Vec<f32> = (0..angles.len()).map(|i| i as f32).collect();

        // polynomial with terms [c0, c1, c2, ..] (c0 + x c1 + x^2 c2)
        let mut polynomial: [f32; 3] = polyfit_rs::polyfit_rs::polyfit(&angles, &indices, 2)?
            .try_into()
            .map_err(|_| "Failed to convert polynomial coefficients")?;

        let errors: Vec<f32> = angles
            .iter()
            .enumerate()
            .map(|x| polyval(&polynomial, *x.1) - (x.0 as f32))
            .collect();

        // let min = errors.iter().cloned().fold(f32::INFINITY, f32::min);
        let max = errors.iter().cloned().fold(f32::NEG_INFINITY, f32::max);

        polynomial[0] -= max;

        let pairs = index_pairs.iter().map(|x| x.0).collect();

        Ok(InterStarIndex {
            pairs,
            angles,
            polynomial,
        })
    }

    /// Get star pairs that match given inter star angle.
    fn pair_lookup(&self, inter_star_angle: f32, tolerance_angle: f32) -> &[[u32; 2]] {
        let lower_threshold = maths_rs::max(inter_star_angle - tolerance_angle, 0.0);
        let upper_threshold = maths_rs::max(inter_star_angle + tolerance_angle, 0.0);
        let lower_index_float = polyval(&self.polynomial, lower_threshold);
        let upper_index_float = polyval(&self.polynomial, upper_threshold);
        let max = self.angles.len() - 1;
        let mut lower_index = (lower_index_float as usize).clamp(0, max);
        let mut upper_index = (upper_index_float as usize).clamp(0, max);

        lower_index = maths_rs::min(lower_index, upper_index);

        while lower_index < max {
            if self.angles[lower_index] > lower_threshold {
                break;
            } else {
                lower_index += 1;
            }
        }
        while upper_index < max {
            if self.angles[upper_index] > upper_threshold {
                break;
            } else {
                upper_index += 1;
            }
        }

        &self.pairs[lower_index..upper_index]
    }
}

#[derive(Debug)]
pub struct MatchResult {
    pub quat: [f32; 4],
    pub match_ids: Vec<u32>,
    pub n_matches: u32,
    pub obs_matched: Vec<[f32; 3]>,
    pub obs_indices: Vec<u32>,
}

#[derive(Debug)]
pub enum FailureReason {
    Unspecified,
    Timeout,
    NotEnoughStars,
    SearchExhausted,
}

#[derive(Debug)]
pub struct DiagnosticData {
    pub reason: FailureReason,
    pub first_svd_failures: usize,
    pub second_svd_failures: usize,
    pub third_svd_failures: usize,
    pub max_first_matches: usize,
    pub max_refined_matches: usize,
}

impl DiagnosticData {
    pub fn new() -> Self {
        DiagnosticData {
            reason: FailureReason::Unspecified,
            first_svd_failures: 0,
            second_svd_failures: 0,
            third_svd_failures: 0,
            max_first_matches: 0,
            max_refined_matches: 0,
        }
    }
}

impl fmt::Display for DiagnosticData {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Delegate to Debug
        write!(f, "{:?}", self)
    }
}

pub struct StarMatcher {
    stars_xyz: Vec<[f32; 3]>,
    star_index: tree::UnitVectorLookup,
    inter_star_index: InterStarIndex,
    /// tolerance of inter star angle in rad
    inter_star_angle_tolerance: f32,
    max_inter_star_angle: f32,
    n_minimum_matches: usize,
    timeout_secs: f32,
}

impl StarMatcher {
    pub fn new(
        stars_xyz: Vec<[f32; 3]>,
        stars_mag: &[f32],
        max_lookup_magnitude: f32,
        max_inter_star_angle: f32,
        inter_star_angle_tolerance: f32,
        n_minimum_matches: usize,
        timeout_secs: f32,
    ) -> Result<Self, &'static str> {
        let star_index = tree::UnitVectorLookup::new(&stars_xyz);
        let inter_star_index = InterStarIndex::new(
            &star_index,
            &stars_xyz,
            &stars_mag,
            max_inter_star_angle,
            max_lookup_magnitude,
        )?;
        Ok(StarMatcher {
            stars_xyz,
            star_index,
            inter_star_index,
            inter_star_angle_tolerance,
            max_inter_star_angle,
            n_minimum_matches,
            timeout_secs,
        })
    }

    pub fn find(&self, obs_xyz: &[[f32; 3]]) -> Result<MatchResult, DiagnosticData> {
        let start_instant = Instant::now();

        let mut diagnostic_data = DiagnosticData::new();

        let iter = match self.triangle_combinations_iterator(obs_xyz.len() as u32) {
            Ok(x) => x,
            Err(_) => {
                diagnostic_data.reason = FailureReason::NotEnoughStars;
                return Err(diagnostic_data);
            }
        };

        for obs_indices in iter {
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
            let ab_pairs = self
                .inter_star_index
                .pair_lookup(angle_ab, self.inter_star_angle_tolerance);
            let ac_pairs = self
                .inter_star_index
                .pair_lookup(angle_ac, self.inter_star_angle_tolerance);
            let bc_pairs = self
                .inter_star_index
                .pair_lookup(angle_bc, self.inter_star_angle_tolerance);

            let iter_finder = trianglefinder::IterTriangleFinder::new(ab_pairs, ac_pairs, bc_pairs);

            // Iterate over possible matching triangles
            for value in iter_finder {
                if let Some(x) = self.check(&obs_indices, &obs_xyz, &value, &mut diagnostic_data) {
                    return Ok(x);
                };
            }

            if start_instant.elapsed().as_secs_f32() > self.timeout_secs {
                diagnostic_data.reason = FailureReason::Timeout;
                return Err(diagnostic_data);
            }
        }
        diagnostic_data.reason = FailureReason::SearchExhausted;
        Err(diagnostic_data)
    }

    /// Iterator over combinations of stars forming triangles.
    fn triangle_combinations_iterator(
        &self,
        n: u32,
    ) -> Result<impl Iterator<Item = [u32; 3]>, &'static str> {
        crate::ordered_combinations::OrderedCombinations::<3>::new(n)
    }

    fn check(
        &self,
        obs_indices: &[u32; 3],
        obs_xyz: &[[f32; 3]],
        cat_indices: &[u32; 3],
        diagnostic_data: &mut DiagnosticData,
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
                diagnostic_data.first_svd_failures += 1;
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
        for obs_i in 0..obs_transformed.ncols() {
            // version 1: use spatial look up
            let obs_vec = obs_transformed.column(obs_i);
            let obs = [obs_vec[0], obs_vec[1], obs_vec[2]];

            // Look up closest star in the catalog to the transformed position of the observation
            let closest_index = self.star_index.lookup_nearest(&obs);

            // Use star if it is close than the allowed threshold
            let closest_cat_star = &self.stars_xyz[closest_index];
            if tree::dot_product(&obs, closest_cat_star) >= dotp_threshold {
                selected_cat_xyz.push(*closest_cat_star);
                selected_obs_xyz.push(obs_xyz[obs_i]);
            }
        }

        // Do not proceed if there are less than the minimum required amount of stars
        if selected_cat_xyz.len() < self.n_minimum_matches {
            diagnostic_data.max_first_matches =
                usize::max(diagnostic_data.max_first_matches, selected_cat_xyz.len());
            return None;
        }

        // Fit rotation matrix on selected observations
        let rotm = match attitude_svd(&selected_cat_xyz, &selected_obs_xyz) {
            None => {
                diagnostic_data.second_svd_failures += 1;
                return None;
            }
            Some(value) => value.cast::<f32>(),
        };

        // Transform all observations, such that they align with the catalog stars
        let obs_xyz_mat =
            na::Matrix3xX::from_iterator(obs_xyz.len(), obs_xyz.into_iter().flatten().cloned());
        let obs_transformed = rotm * obs_xyz_mat;

        // Do another pass with the improved rotation matrix
        selected_obs_xyz.clear();
        selected_cat_xyz.clear();
        let mut selected_cat_indices = Vec::new();
        let mut selected_obs_indices = Vec::new();
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
                selected_cat_indices.push(closest_index as u32);
                selected_obs_indices.push(obs_i as u32);
            }
        }

        // Do not proceed if there are less than the minimum required amount of stars
        if selected_cat_xyz.len() < self.n_minimum_matches {
            diagnostic_data.max_refined_matches =
                usize::max(diagnostic_data.max_refined_matches, selected_cat_xyz.len());
            return None;
        }

        // Fit rotation matrix on selected observations
        let final_rotm = match attitude_svd(&selected_cat_xyz, &selected_obs_xyz) {
            None => {
                diagnostic_data.third_svd_failures += 1;
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
            obs_matched: selected_obs_xyz,
            obs_indices: selected_obs_indices,
        })
    }

    pub fn stars_xyz(&self) -> &[[f32; 3]] {
        &self.stars_xyz
    }
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
            (a * d) as f64,
            (b * d) as f64,
            (c * d) as f64,
            (a * e) as f64,
            (b * e) as f64,
            (c * e) as f64,
            (a * f) as f64,
            (b * f) as f64,
            (c * f) as f64,
        );
        mat += outer_prod;
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
