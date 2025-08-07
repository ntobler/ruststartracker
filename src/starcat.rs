use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use chrono::{DateTime, TimeZone, Utc};
use csv::ReaderBuilder;
use nalgebra::Vector3;
use serde::Deserialize;

const AU: f64 = 149_597_870.693;
const GAIA_EPOCH: f64 = 2016.0;
#[cfg(feature = "gaia")]
const GAIA_2016_CSV: &str = include_str!("../ruststartracker/gaia_data_j2016.csv");

pub trait Cast<T> {
    fn cast(self) -> T;
}

impl Cast<f32> for f64 {
    #[inline(always)]
    fn cast(self) -> f32 {
        self as f32
    }
}

impl Cast<f64> for f64 {
    #[inline(always)]
    fn cast(self) -> f64 {
        self
    }
}

#[derive(Debug, Deserialize)]
pub struct StarRaw {
    pub source_id: u64,
    pub ra: f64,
    pub dec: f64,
    pub parallax: f64,
    pub pmra: f64,
    pub pmdec: f64,
    pub phot_g_mean_mag: f64,
}

pub struct Star {
    pub source_id: u64,
    pub ra: f64,
    pub dec: f64,
    pub parallax: f64,
    pub proper_motion_ra: f64,
    pub proper_motion_dec: f64,
    pub magnitude: f64,
}

impl Star {
    const DEG2RAD: f64 = core::f64::consts::PI / 180.0;
    const MAS2RAD: f64 = core::f64::consts::PI / (180.0 * 3600.0 * 1000.0);

    pub fn from_raw(raw: StarRaw) -> Self {
        Star {
            source_id: raw.source_id,
            ra: raw.ra * Star::DEG2RAD,
            dec: raw.dec * Star::DEG2RAD,
            parallax: raw.parallax * Star::MAS2RAD,
            proper_motion_ra: raw.pmra * Star::MAS2RAD,
            proper_motion_dec: raw.pmdec * Star::MAS2RAD,
            magnitude: raw.phot_g_mean_mag,
        }
    }
}

fn time_to_epoch(t: DateTime<Utc>) -> f64 {
    let delta_t_j2000 = 64.0;
    let j2000 = Utc
        .with_ymd_and_hms(2000, 1, 1, 12, 0, 0)
        .unwrap()
        .timestamp() as f64
        - delta_t_j2000;
    let seconds_in_j_year = 365.25 * 86400.0;
    (t.timestamp() as f64 - j2000) / seconds_in_j_year + 2000.0
}

pub struct StarCatalog {
    pub stars: Vec<Star>,
    pub epoch: f64,
}

impl StarCatalog {
    pub fn new_from_file<P: AsRef<Path>>(
        filename: P,
        epoch: f64,
        max_magnitude: Option<f64>,
    ) -> Result<Self, String> {
        let file = File::open(filename).map_err(|_| "Unable to open star catalog file")?;
        let reader = BufReader::new(file);
        Self::new_from_buffer(reader, epoch, max_magnitude)
    }

    pub fn new_from_string(
        data: &str,
        epoch: f64,
        max_magnitude: Option<f64>,
    ) -> Result<Self, String> {
        let reader = BufReader::new(data.as_bytes());
        Self::new_from_buffer(reader, epoch, max_magnitude)
    }

    #[cfg(feature = "gaia")]
    pub fn new_from_gaia(max_magnitude: Option<f64>) -> Result<Self, String> {
        StarCatalog::new_from_string(GAIA_2016_CSV, GAIA_EPOCH, max_magnitude)
    }

    fn new_from_buffer<T: std::io::Read>(
        reader: BufReader<T>,
        epoch: f64,
        max_magnitude: Option<f64>,
    ) -> Result<Self, String> {
        let mut csv_reader = ReaderBuilder::new().from_reader(reader);
        let mut stars = Vec::new();
        for line in csv_reader.deserialize() {
            let star: Star = Star::from_raw(line.map_err(|e| e.to_string())?);
            if star.magnitude < max_magnitude.unwrap_or(f64::INFINITY) {
                stars.push(star);
            }
        }
        Ok(StarCatalog { stars, epoch })
    }

    pub fn normalized_positions<T>(
        &self,
        epoch: Option<f64>,
        observer_position: Option<[f64; 3]>,
    ) -> Vec<[T; 3]>
    where
        f64: Cast<T>,
    {
        let epoch = epoch.unwrap_or_else(|| time_to_epoch(Utc::now()));
        let delta_epoch = epoch - self.epoch;

        let mut vectors: Vec<[T; 3]> = Vec::with_capacity(self.stars.len());

        let parallax_correction_factor = match observer_position {
            Some(obs_pos) => {
                let obs_pos_vec = Vector3::from_column_slice(&obs_pos);
                Some(obs_pos_vec / AU)
            }
            None => None,
        };

        for star in self.stars.iter() {
            let cos_ra = star.ra.cos();
            let sin_ra = star.ra.sin();
            let cos_dec = star.dec.cos();
            let sin_dec = star.dec.sin();

            let mut vec = Vector3::new(cos_dec * cos_ra, cos_dec * sin_ra, sin_dec);

            // Proper motion correction
            let p_hat = Vector3::new(-sin_ra, cos_ra, 0.0);
            let q_hat = Vector3::new(-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec);
            let pm = delta_epoch * (star.proper_motion_ra * p_hat + star.proper_motion_dec * q_hat);
            vec += pm;

            // Parallax correction
            if let Some(f) = parallax_correction_factor {
                let plx = star.parallax * f;
                vec -= plx;
            }

            // Normalize
            let norm_vec: [f64; 3] = vec.normalize().into();

            vectors.push([
                Cast::<T>::cast(norm_vec[0]),
                Cast::<T>::cast(norm_vec[1]),
                Cast::<T>::cast(norm_vec[2]),
            ]);
        }

        vectors
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_starcat() {
        let cat = StarCatalog::new_from_file(
            "ruststartracker/gaia_data_j2016.csv",
            GAIA_EPOCH,
            Some(5.5),
        )
        .unwrap();

        assert!(cat.epoch == GAIA_EPOCH);
        assert!(cat.stars.len() > 1000);

        let positions: Vec<[f64; 3]> = cat.normalized_positions(Some(2025.23), None);
        println!("{:?}", positions.len());
        let positions: Vec<[f32; 3]> = cat.normalized_positions(Some(2025.23), None);
        println!("{:?}", positions.len());
    }

    #[test]
    fn test_proper_motion_star() {
        let cat = StarCatalog::new_from_file(
            "ruststartracker/gaia_data_j2016.csv",
            GAIA_EPOCH,
            Some(8.0),
        )
        .unwrap();

        let id = 725422076533653504;

        let index = cat.stars.iter().position(|s| s.source_id == id).unwrap();

        let pos_2010: [f64; 3] = cat.normalized_positions(Some(2010.0), None)[index];
        let pos_2030: [f64; 3] = cat.normalized_positions(Some(2030.0), None)[index];

        let angle_rad = f64::acos(
            pos_2010[0] * pos_2030[0] + pos_2010[1] * pos_2030[1] + pos_2010[2] * pos_2030[2],
        );

        let angle_mas = angle_rad / Star::MAS2RAD;

        let expected_travel_mas = 20.0 * 22.717;

        assert!(
            (angle_mas / expected_travel_mas).abs() - 1.0 < 0.001,
            "Expected travel: {}, got: {}",
            expected_travel_mas,
            angle_mas
        );
    }

    #[cfg(feature = "gaia")]
    #[test]
    fn test_gaia() {
        let cat1 = StarCatalog::new_from_file(
            "ruststartracker/gaia_data_j2016.csv",
            GAIA_EPOCH,
            Some(5.5),
        )
        .unwrap();
        let cat2 = StarCatalog::new_from_gaia(Some(5.5)).unwrap();
        assert!(cat1.stars.len() == cat2.stars.len());
    }
}
