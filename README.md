# Lightweight Python Star Tracker With Rust Backend

Based on the methodology used in https://github.com/nasa/COTS-Star-Tracker, with following improvements:
- Reduced dependencies to opencv and numpy for lightweight usage in a Raspberry Pi.
- Reimplemented computationally expensive parts in rust. This includes most parts that are not image processing related.
- Added quadratic inter star angle index look up polynomial for faster triangle search.
- Added spatial index to look up neighboring stars.

Features:
- Attitude estimation from image and camera calibration parameters.
- Attitude estimation from list of star observation coordinates.
- Star catalog creation with temporal corrections.

## Example

### Rust

See [examples/basic.rs](examples/basic.rs)

```rust
// Get catalog positions
let catalog: StarCatalog = StarCatalog::from_gaia(max_magnitude: ...).unwrap();
let stars_xyz: Vec<[f32; 3]> = catalog.normalized_positions(epoch: ..., observer_position: ...);
let stars_mag: Vec<f32> = catalog.magnitudes();

// Create StarTracker instance (reuse this)
let star_matcher = StarMatcher::new(
    stars_xyz,
    stars_mag,
    max_lookup_magnitude: ...
    max_inter_star_angle: ...,
    inter_star_angle_tolerance: ...,
    min_matches: ...,
    timeout: ...
);

// Normalized observation in the camera frame
let obs_xyz_camera: Vec<[f32; 3]> = ...

let result = star_matcher.find(&obs_xyz_camera);
println!("Result: {:?}", result);
```

### Python
```python
import ruststartracker

# Get catalog positions
catalog = ruststartracker.StarCatalog(max_magnitude=...)
star_catalog_vecs = catalog.normalized_positions(epoch=...)

# Define opencv camera parameters, see https://docs.opencv.org/4.x/dc/dbb/tutorial_py_calibration.html
camera_params = ruststartracker.CameraParameters(
    camera_matrix=...,
    cam_resolution=...,
    dist_coefs=...,
)

# Create StarTracker instance (reuse this)
st = ruststartracker.StarTracker(
    star_catalog_vecs,
    camera_params,
    max_inter_star_angle=...,
    inter_star_angle_tolerance=...,
    n_minimum_matches=...,
    timeout_secs=...,
)

# Obtain numpy array image
img = ...

# Find attitude from given image
result = st.process_image(img)

print(result)
# StarTrackerResult(quat=[-0.43977802991867065, -0.439766526222229, -0.4398997128009796, 0.6478340029716492], match_ids=[1435, 1272, 1140, 2035, 1070, 1438, 1338, 903, 260, 2141, 1771, 1727, 385, 1717, 2204, 2062, 1989, 1634, 708, 1357], n_matches=20, duration_s=0.0003700880042742938)
```

## Installation

- Install with `pip install ruststartracker` (Currently only ARM/x86 Linux wheels available).

## Attributions

### Gaia Data

This project includes data from the European Space Agency (ESA) mission [**Gaia**](https://www.cosmos.esa.int/gaia), processed by the **Gaia Data Processing and Analysis Consortium (DPAC)**.
Funding for the DPAC has been provided by national institutions, in particular the institutions participating in the Gaia Multilateral Agreement.

Gaia DR3 data is © European Space Agency and is released under the [**Creative Commons Attribution 4.0 International License (CC BY 4.0)**](https://creativecommons.org/licenses/by/4.0/).

> Gaia Collaboration, Vallenari et al. (2022), *A\&A* **674**, A1.
> [DOI: 10.1051/0004-6361/202243940](https://doi.org/10.1051/0004-6361/202243940)
