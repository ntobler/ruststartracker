use std::fmt;

use crate::tree;
use crate::trianglefinder;
use std::time::Instant;

extern crate nalgebra as na;

/// Naive implementation for tree::UnitVectorLookup::look_up_close_angles
/// Calculates all star angles, instead of looking up neighbors in a spatial index
/// Turned out to be faster in some cases
pub fn look_up_close_angles_naive(
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
        for b in (a + 1)..vectors.len() {
            if magnitudes[b] > max_magnitude {
                continue;
            }
            let vec_b = &vectors[b];
            let dotp = tree::dot_product(vec_a, vec_b);
            if dotp >= cos_max_angle {
                index_pairs.push(([a as u32, b as u32], dotp));
            }
        }
    }
    index_pairs
}

pub struct InterStarIndex {
    pub pairs: Vec<[u32; 2]>,
    pub cos_angles: Vec<f32>,
    // polynomial with terms [c0, c1, c2] (c0 + x c1 + x^2 c2)
    pub polynomial: [f32; 3],
}

impl InterStarIndex {
    pub fn new(
        star_index: &tree::UnitVectorLookup,
        stars_xyz: &[[f32; 3]],
        stars_mag: &[f32],
        cos_max_angle: f32,
        max_magnitude: f32,
    ) -> Result<InterStarIndex, &'static str> {
        if stars_xyz.len() != stars_mag.len() {
            return Err("stars_xyz and stars_mag must have same length");
        }

        let mut index_pairs =
            star_index.look_up_close_angles(stars_xyz, stars_mag, cos_max_angle, max_magnitude);

        if index_pairs.is_empty() {
            return Err("Given star positions do not result in any angles below the threshold.");
        }

        // Sort pairs by increasing angle (cosine of angle is decreasing)
        index_pairs.sort_unstable_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

        let scale = (index_pairs.len() + 1) as f32;
        let inv_scale = 1.0 / scale;
        let scaled_indices: Vec<f32> = (0..index_pairs.len())
            .map(|i| i as f32 * inv_scale)
            .collect();
        let pairs: Vec<[u32; 2]> = index_pairs.iter().map(|x| x.0).collect();
        let cos_angles: Vec<f32> = index_pairs.iter().map(|x| x.1).collect();

        let transformed_cos_angles: Vec<f32> = cos_angles
            .iter()
            .map(|cos_angle: &f32| Self::transform_angle(*cos_angle))
            .collect();

        // Fit polynomial
        // polynomial with terms [c0, c1, c2, ..] (c0 + x c1 + x^2 c2)
        let mut polynomial: [f32; 3] =
            polyfit_rs::polyfit_rs::polyfit(&transformed_cos_angles, &scaled_indices, 2)?
                .try_into()
                .map_err(|_| "Failed to convert polynomial coefficients")?;

        let diff_poly_to_value: Vec<f32> = transformed_cos_angles
            .iter()
            .enumerate()
            .map(|(i, &x)| Self::polyval(&polynomial, x) * scale - (i as f32))
            .collect();

        // let min = diff_poly_to_value
        //     .iter()
        //     .cloned()
        //     .fold(f32::INFINITY, f32::min);
        let max = diff_poly_to_value
            .iter()
            .cloned()
            .fold(f32::NEG_INFINITY, f32::max);

        polynomial[0] -= max * inv_scale;

        Ok(InterStarIndex {
            pairs,
            cos_angles,
            polynomial,
        })
    }

    /// Transform cosine of angle to a value that is more suitable for polynomial fitting
    fn transform_angle(cos_angle: f32) -> f32 {
        // 10 is an arbitrary scaling factor that maps an angle of approx 25 degrees to a value of 1.0
        (1.0 - cos_angle) * 10.0
    }

    /// Evaluate a polynomial
    fn polyval<const S: usize>(coeffs: &[f32; S], x: f32) -> f32 {
        coeffs.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    /// Get star pairs that match given inter star angle.
    pub fn pair_lookup(&self, cos_inter_star_angle: f32, cos_tolerance_angle: f32) -> &[[u32; 2]] {
        // Both sine and cosine of the inter star angle are positive if it is between 0 and 1
        let c = cos_inter_star_angle;
        let s = (1.0 - c * c).sqrt();

        // Both sine and cosine of the tolerance are positive if it is between 0 and 1
        let tc = cos_tolerance_angle;
        let ts = (1.0 - tc * tc).sqrt();

        // Essentially equal to cos(acos(c) +/- acos(tc)), but faster to compute
        let cos_lower_threshold = c * tc + s * ts;
        let cos_upper_threshold = c * tc - s * ts;

        let scale = (self.pairs.len() + 1) as f32;
        let lower_index_float =
            Self::polyval(&self.polynomial, Self::transform_angle(cos_lower_threshold)) * scale;
        let upper_index_float =
            Self::polyval(&self.polynomial, Self::transform_angle(cos_upper_threshold)) * scale;
        let max = self.cos_angles.len() - 1;
        let mut lower_index = (lower_index_float as usize).clamp(0, max);
        let mut upper_index = (upper_index_float as usize).clamp(0, max);

        lower_index = usize::min(lower_index, upper_index);

        while lower_index < max {
            if self.cos_angles[lower_index] < cos_lower_threshold {
                break;
            } else {
                lower_index += 1;
            }
        }
        while upper_index < max {
            if self.cos_angles[upper_index] < cos_upper_threshold {
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
    pub svd_failures: usize,
    pub max_first_matches: usize,
    pub max_refined_matches: usize,
}

impl DiagnosticData {
    pub fn new() -> Self {
        DiagnosticData {
            reason: FailureReason::Unspecified,
            svd_failures: 0,
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
    cos_inter_star_angle_tolerance: f32,
    max_cos_inter_star_angle: f32,
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
            max_inter_star_angle.cos(),
            max_lookup_magnitude,
        )?;
        Ok(StarMatcher {
            stars_xyz,
            star_index,
            inter_star_index,
            cos_inter_star_angle_tolerance: inter_star_angle_tolerance.cos(),
            max_cos_inter_star_angle: max_inter_star_angle.cos(),
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
            let cos_angle_ab = tree::dot_product(&obs_a, &obs_b);
            let cos_angle_ac = tree::dot_product(&obs_a, &obs_c);
            let cos_angle_bc = tree::dot_product(&obs_b, &obs_c);

            // Continue if any angle is larger than the maximum
            // angle contained in the inter star angle index
            let t = self.max_cos_inter_star_angle;
            if (cos_angle_ab < t) || (cos_angle_ac < t) || (cos_angle_bc < t) {
                continue;
            }

            // Look up all pairs in the catalog that could match the observations
            let ab_pairs = self
                .inter_star_index
                .pair_lookup(cos_angle_ab, self.cos_inter_star_angle_tolerance);
            let ac_pairs = self
                .inter_star_index
                .pair_lookup(cos_angle_ac, self.cos_inter_star_angle_tolerance);
            let bc_pairs = self
                .inter_star_index
                .pair_lookup(cos_angle_bc, self.cos_inter_star_angle_tolerance);

            let iter_finder = trianglefinder::IterTriangleFinder::new(ab_pairs, ac_pairs, bc_pairs);

            // Iterate over possible matching triangles
            for abc_candidate in iter_finder {
                match self.check(&obs_indices, &obs_xyz, &abc_candidate) {
                    Ok(x) => return Ok(x),
                    Err(StarError::SVDFailure) => {
                        diagnostic_data.svd_failures += 1;
                    }
                    Err(StarError::NotEnoughFirstMatches(n)) => {
                        diagnostic_data.max_first_matches =
                            usize::max(diagnostic_data.max_first_matches, n);
                    }
                    Err(StarError::NotEnoughRefinedMatches(n)) => {
                        diagnostic_data.max_refined_matches =
                            usize::max(diagnostic_data.max_refined_matches, n);
                    }
                    Err(_) => {
                        // Other errors are not tracked in the diagnostics
                    }
                }
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
    ) -> Result<MatchResult, StarError> {
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
        let rotm = attitude_svd(&cat_triangle_xyz, &obs_triangle_xyz)?.cast::<f32>();

        // Transform all observations, such that they align with the catalog stars
        let obs_xyz_mat =
            na::Matrix3xX::from_iterator(obs_xyz.len(), obs_xyz.into_iter().flatten().cloned());
        let obs_transformed = rotm * obs_xyz_mat;

        // Find close neighbors for each observation
        let dotp_threshold = self.cos_inter_star_angle_tolerance;
        let mut selected_obs_xyz = Vec::new();
        let mut selected_cat_xyz = Vec::new();
        for obs_i in 0..obs_transformed.ncols() {
            // version 1: use spatial look up
            let obs_vec = obs_transformed.column(obs_i);
            let obs = std::array::from_fn(|i| obs_vec[i]);

            // Look up closest star in the catalog to the transformed position of the observation
            let closest_index = self.star_index.lookup_nearest(&obs);

            // Use star if it is closer than the allowed threshold
            let closest_cat_star = &self.stars_xyz[closest_index];
            if tree::dot_product(&obs, closest_cat_star) >= dotp_threshold {
                selected_cat_xyz.push(*closest_cat_star);
                selected_obs_xyz.push(obs_xyz[obs_i]);
            }
        }

        // Do not proceed if there are less than the minimum required amount of stars
        if selected_cat_xyz.len() < self.n_minimum_matches {
            return Err(StarError::NotEnoughFirstMatches(selected_cat_xyz.len()));
        }

        // Fit rotation matrix on selected observations
        let rotm = attitude_svd(&selected_cat_xyz, &selected_obs_xyz)?.cast::<f32>();

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
            return Err(StarError::NotEnoughRefinedMatches(selected_cat_xyz.len()));
        }

        // Fit rotation matrix on selected observations
        let final_rotm = attitude_svd(&selected_cat_xyz, &selected_obs_xyz)?;

        // Convert rotation matrix into quaternion
        let quat = na::UnitQuaternion::from_rotation_matrix(&na::Rotation3::from_matrix_unchecked(
            final_rotm,
        ));

        // Return result with statistics
        Ok(MatchResult {
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

#[derive(thiserror::Error, Debug)]
pub enum StarError {
    #[error("Failed to create inter star index: {0}")]
    InterStarIndexCreationError(String),
    #[error("Failed to calculate singular value decomposition")]
    SVDFailure,
    #[error("Not enough first matches {0}")]
    NotEnoughFirstMatches(usize),
    #[error("Not enough refined matches {0}")]
    NotEnoughRefinedMatches(usize),
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
pub fn attitude_svd(
    cat_xyz: &[[f32; 3]],
    obs_xyz: &[[f32; 3]],
) -> Result<na::Matrix3<f64>, StarError> {
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
        true => Ok(t),
        false => Err(StarError::SVDFailure),
    }
}
