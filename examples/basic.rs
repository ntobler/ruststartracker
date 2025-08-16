use ruststartracker::star::StarMatcher;
use ruststartracker::starcat::StarCatalog;

fn get_observations() -> Vec<[f32; 3]> {
    // For demo purposes we extract some bright stars from the catalog.

    // Get catalog positions
    let catalog = StarCatalog::new_from_gaia(Some(5.0)).unwrap();
    let stars_xyz: Vec<[f32; 3]> = catalog.normalized_positions(Some(2025.0), None);

    // Get some observations from the catalog
    stars_xyz
        .iter()
        .filter(|x| x[1] > f32::cos(0.5))
        .map(|x| *x)
        .collect()
}

fn main() {
    // Get catalog positions
    let catalog = StarCatalog::new_from_gaia(Some(6.0)).unwrap();
    let stars_xyz: Vec<[f32; 3]> = catalog.normalized_positions(Some(2025.0), None);
    let stars_mag: Vec<f32> = catalog.magnitudes();

    // Create StarTracker instance (reuse this)
    let star_matcher = StarMatcher::new(stars_xyz, &stars_mag, 5.0, 1.0, 0.002, 10, 0.2).unwrap();

    // Get observation in the camera frame (provide this function)
    let obs_xyz_camera: Vec<[f32; 3]> = get_observations();

    // Lookup attitude
    let result = star_matcher.find(&obs_xyz_camera).unwrap();

    println!("Result: {:?}", result);
}
