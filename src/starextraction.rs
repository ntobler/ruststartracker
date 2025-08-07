#[cfg(feature = "improc")]
use opencv;
use opencv::{core::Mat, prelude::*};

pub fn get_threshold_from_histogram(image_row_major: &[u8], fraction: f64) -> u8 {
    // Calculate histogram
    let mut hist: [usize; 256] = [0; 256];
    for &v in image_row_major {
        hist[v as usize] += 1;
    }

    // Calculate cumulative sum of histogram
    let mut acc = 0;
    for i in 0..hist.len() {
        acc += hist[i];
        hist[i] = acc;
    }

    // Define threshold in terms of cumulative sum value
    let total: usize = *hist.last().unwrap(); // Safe: hist always has 256 elements
    let threshold_value = (total as f64 * f64::clamp(fraction, 0.0, 1.0)) as usize;

    //Find index in histogram where threshold is exceeded
    // Find threshold index
    hist.iter()
        .position(|&v| v >= threshold_value)
        .unwrap_or(255) as u8
}

pub fn threshold(image: &[u8], threshold: u8) -> Vec<u8> {
    image.iter().map(|&x| (x > threshold) as u8).collect()
}

pub fn extract_observations(
    image_row_major: &[u8],
    imsize: (usize, usize),
    threshold_value: u8,
    min_area: usize,
    max_area: usize,
) -> Result<(Vec<[f64; 2]>, Vec<f64>), String> {
    let (width, height) = imsize;

    let thresholded_row_major: Vec<u8> = threshold(image_row_major, threshold_value);

    // Call opencv function
    let thresholded = opencv::core::Mat::new_rows_cols_with_data(
        height as i32,
        width as i32,
        &thresholded_row_major,
    )
    .map_err(|e| format!("Failed to create Mat: {}", e))?;
    let mut labels = Mat::default();
    let mut stats = Mat::default();
    let mut centroids = Mat::default();
    _ = opencv::imgproc::connected_components_with_stats(
        &thresholded,
        &mut labels,
        &mut stats,
        &mut centroids,
        4,
        opencv::core::CV_16U,
    )
    .map_err(|e| format!("Failed to perform connected component analysis: {}", e))?;

    let stats_data: &[i32] = stats
        .data_typed::<i32>()
        .map_err(|e| format!("Failed to get stats data: {}", e))?;
    let stats_rows = stats.rows() as usize;
    let stats_cols = stats.cols() as usize;

    let mut centers_intensities = Vec::new();

    for i in 0..stats_rows {
        // Get bounding box around star
        let stats_base = i * stats_cols;
        let x = stats_data[stats_base + opencv::imgproc::CC_STAT_LEFT as usize] as usize;
        let y = stats_data[stats_base + opencv::imgproc::CC_STAT_TOP as usize] as usize;
        let w = stats_data[stats_base + opencv::imgproc::CC_STAT_WIDTH as usize] as usize;
        let h = stats_data[stats_base + opencv::imgproc::CC_STAT_HEIGHT as usize] as usize;
        let area = stats_data[stats_base + opencv::imgproc::CC_STAT_AREA as usize] as usize;

        if area < min_area
            || area > max_area
            || x == 0
            || y == 0
            || x + w >= width
            || y + h >= height
        {
            continue;
        }

        let mut mid_x = 0;
        let mut mid_y = 0;
        let mut intensity_sum = 0;

        // Accumulate intensity and intensity weighted positions in the window
        for yy in y - 1..y + h + 1 {
            //row should be inner loop for better cache optimization
            let img_base = yy * width;
            for xx in x - 1..x + w + 1 {
                let value = image_row_major[img_base + xx] as usize;
                intensity_sum += value;
                mid_x += xx * value;
                mid_y += yy * value;
            }
        }

        if intensity_sum == 0 {
            continue;
        }

        centers_intensities.push((
            mid_x as f64 / intensity_sum as f64,
            mid_y as f64 / intensity_sum as f64,
            intensity_sum as f64,
        ));
    }

    // Sort centers by intensity in descending order
    centers_intensities.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
    let centers = centers_intensities.iter().map(|x| [x.0, x.1]).collect();
    let intensities = centers_intensities.iter().map(|x| x.2).collect();

    Ok((centers, intensities))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple() {
        #[rustfmt::skip]
        let data= [
            0, 0, 0, 0, 0, 0, 0,
            0, 4, 4, 0, 0, 0, 0,
            0, 4, 4, 0, 4, 8, 0,
            0, 0, 0, 0, 4, 8, 0,
            0, 0, 0, 0, 0, 0, 0,
        ];

        let (centers, intensities) = extract_observations(&data, (7, 5), 3, 2, 6).unwrap();

        let centers_expected = vec![[4.666666, 2.5], [1.5, 1.5]];
        let intensities_expected = vec![24.0, 16.0];

        // Helper function to compare floats with tolerance
        fn assert_vec2_close(a: &[f64; 2], b: &[f64; 2], tol: f64) {
            assert!((a[0] - b[0]).abs() < tol, "x: {} vs {}", a[0], b[0]);
            assert!((a[1] - b[1]).abs() < tol, "y: {} vs {}", a[1], b[1]);
        }

        fn assert_vec_close(a: &[f64], b: &[f64], tol: f64) {
            assert_eq!(a.len(), b.len());
            for (x, y) in a.iter().zip(b.iter()) {
                assert!((x - y).abs() < tol, "{} vs {}", x, y);
            }
        }

        assert_eq!(centers.len(), centers_expected.len());
        for (c, ce) in centers.iter().zip(centers_expected.iter()) {
            assert_vec2_close(c, ce, 1e-4);
        }
        assert_vec_close(&intensities, &intensities_expected, 1e-4);

        // Additional test: single bright pixel
        #[rustfmt::skip]
        let data2 = [
            0, 0, 0,
            0, 9, 0,
            0, 0, 0
        ];
        let (centers2, intensities2) = extract_observations(&data2, (3, 3), 5, 1, 2).unwrap();
        let centers2_expected = vec![[1.0, 1.0]];
        let intensities2_expected = vec![9.0];
        assert_eq!(centers2.len(), centers2_expected.len());
        for (c, ce) in centers2.iter().zip(centers2_expected.iter()) {
            assert_vec2_close(c, ce, 1e-4);
        }
        assert_vec_close(&intensities2, &intensities2_expected, 1e-4);

        // Additional test: two separated stars
        #[rustfmt::skip]
        let data3 = [
            0, 0, 0, 0, 0,
            0, 7, 0, 8, 0,
            0, 0, 0, 0, 0
        ];
        let (centers3, intensities3) = extract_observations(&data3, (5, 3), 5, 1, 2).unwrap();
        let centers3_expected = vec![[1.0, 1.0], [3.0, 1.0]];
        let intensities3_expected = vec![7.0, 8.0];
        // Order is not guaranteed, so sort by x
        let mut zipped: Vec<_> = centers3.iter().zip(intensities3.iter()).collect();
        zipped.sort_by(|a, b| a.0[0].partial_cmp(&b.0[0]).unwrap());
        let (centers3, intensities3): (Vec<_>, Vec<_>) =
            zipped.into_iter().map(|(c, i)| (*c, *i)).unzip();
        for (c, ce) in centers3.iter().zip(centers3_expected.iter()) {
            assert_vec2_close(c, ce, 1e-4);
        }
        assert_vec_close(&intensities3, &intensities3_expected, 1e-4);

        // Additional test: no stars above threshold
        #[rustfmt::skip]
        let data4 = [
            0, 0, 0, 0, 0,
            0, 1, 1, 1, 0,
            0, 1, 1, 1, 0,
            0, 1, 1, 1, 0,
            0, 0, 0, 0, 0,
        ];
        let (centers4, intensities4) = extract_observations(&data4, (5, 5), 10, 1, 10).unwrap();
        assert!(centers4.is_empty());
        assert!(intensities4.is_empty());

        // Additional test: stars in all corners
        #[rustfmt::skip]
        let data5 = [
            5, 5, 0, 5, 5,
            5, 5, 1, 5, 5,
            0, 1, 8, 1, 0,
            4, 5, 1, 5, 7,
            4, 4, 0, 7, 7,
        ];
        let (centers5, intensities5) = extract_observations(&data5, (5, 5), 2, 1, 10).unwrap();
        let centers5_expected = vec![[2.0, 2.0]];
        let intensities5_expected = vec![32.0];
        assert_eq!(centers.len(), centers_expected.len());
        for (c, ce) in centers5.iter().zip(centers5_expected.iter()) {
            assert_vec2_close(c, ce, 1e-4);
        }
        assert_vec_close(&intensities5, &intensities5_expected, 1e-4);
    }

    pub fn threshold_opencv(image: &[u8], threshold: u8) -> Mat {
        let image_mat =
            opencv::core::Mat::new_rows_cols_with_data(image.len() as i32, 1, &image).unwrap();
        let mut out = Mat::default();
        opencv::imgproc::threshold(
            &image_mat,
            &mut out,
            threshold as f64,
            1 as f64,
            opencv::imgproc::THRESH_BINARY,
        )
        .unwrap();
        out
    }

    #[test]
    fn test_threshold() {
        use std::time::Instant;

        let image = vec![56; 4000000];

        let start1 = Instant::now();
        let out1 = threshold(&image, 42);
        let duration1 = start1.elapsed();

        let start2 = Instant::now();
        let out2 = threshold_opencv(&image, 42);
        let duration2 = start2.elapsed();

        println!("threshold (Rust) took: {:?}", duration1);
        println!("threshold_opencv (OpenCV) took: {:?}", duration2);

        let sum1: usize = out1.into_iter().map(|x| x as usize).sum();
        let d = out2.data_typed::<u8>().unwrap();
        let sum2: usize = d.into_iter().map(|&x| x as usize).sum();

        println!("{:?}", sum1);
        println!("{:?}", sum2);
    }

    #[test]
    fn test_threshold_vlues() {
        let image = vec![1, 2, 3, 4, 5, 6, 7, 8];

        let out1 = threshold(&image, 4);
        let out2 = threshold_opencv(&image, 4);

        let sum1: usize = out1.into_iter().map(|x| x as usize).sum();
        let d = out2.data_typed::<u8>().unwrap();
        let sum2: usize = d.into_iter().map(|&x| x as usize).sum();

        assert_eq!(sum1, 4);
        assert_eq!(sum2, 4);
    }

    #[test]
    fn test_histogram() {
        let image = vec![0, 4, 4];
        let threshold = get_threshold_from_histogram(&image, 0.9);
        println! {"{:?}", threshold};
        assert_eq!(threshold, 4);

        let image = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9];
        let threshold = get_threshold_from_histogram(&image, 0.899);
        println! {"{:?}", threshold};
        assert_eq!(threshold, 7);
        let threshold = get_threshold_from_histogram(&image, 0.900);
        println! {"{:?}", threshold};
        assert_eq!(threshold, 8);
        let threshold = get_threshold_from_histogram(&image, 0.901);
        println! {"{:?}", threshold};
        assert_eq!(threshold, 8);
    }
}
