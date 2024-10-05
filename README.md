# Leightweight Python Start Tracker With Rust Backend

Based on the methodology used in https://github.com/nasa/COTS-Star-Tracker, with following improvements:
- Reduced dependencies to opencv and numpy for leightweight usage in a Raspberry Pi.
- Reimplemented computationally expensive parts in rust. This includes most parts that are not image processing related.
- Added quadratic inter star angle index look up polynom for faster triangle search.
- Added spatial index to look up neighboring stars.

Features:
- Attitude estimation from image and camera calibration parameters.
- Attitude estimation from list of star obvervation coordinates.
- Star catalog creation with tempral corrections.

## Example

```python
import ruststartracker

# Get catalog positions
catalog = ruststartracker.StarCatalog()
star_catalog_vecs = catalog.normalized_positions()

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
)

# Obtain numpy array image
img = ...

# Find attitude from given image
result = st.process_image(img)

print(result)
# StarTrackerResult(quat=[-0.43977802991867065, -0.439766526222229, -0.4398997128009796, 0.6478340029716492], match_ids=[1435, 1272, 1140, 2035, 1070, 1438, 1338, 903, 260, 2141, 1771, 1727, 385, 1717, 2204, 2062, 1989, 1634, 708, 1357], n_matches=20, duration_s=0.0003700880042742938)
```

## Installation

- Make sure rust tool chain (`cargo`) is installed and in the `PATH` environment variable.
- Install with `pip install git+https://github.com/ntobler/startracker.git`.

## TODOs

- Improve error messages.
- Return more diagnostic data.
