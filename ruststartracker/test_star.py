import time

import numpy as np
import pytest
import scipy.spatial

import ruststartracker
import ruststartracker.libruststartracker
import ruststartracker.star


@pytest.mark.parametrize("impl", ["python", "rust"])
def test_extract_observations(impl: str):
    size_x, size_y = (960, 480)
    img = np.zeros((size_y, size_x), np.uint8)
    points = np.array([(3, 5), (23, 13), (30, 50), (230, 130)])
    for x, y in points:
        img[y - 1 : y + 3, x - 1 : x + 3] = 50

    t0 = time.monotonic()
    if impl == "python":
        centers, intensities = ruststartracker.star._extract_observations(img, threshold=30)
    elif impl == "rust":
        centers, intensities = ruststartracker.libruststartracker.extract_observations(
            img, 30, 3, 300
        )
    else:
        raise AssertionError
    print(f"Extracting observations took {time.monotonic() - t0:.5f} seconds")

    assert isinstance(centers, np.ndarray)
    assert isinstance(intensities, np.ndarray)
    np.testing.assert_almost_equal(centers, points + 0.5)
    np.testing.assert_almost_equal(intensities, 50 * 16)


@pytest.fixture
def setup():
    rng = np.random.default_rng(42)

    n_cat_stars = 2617

    vec = rng.normal(size=[n_cat_stars, 3]).astype(np.float32)
    vec /= np.linalg.norm(vec, axis=-1, keepdims=True)

    mag = rng.uniform(0, 10, size=vec.shape[:1]).astype(np.float32)

    angle_threshold = np.radians(10)
    dotp = np.sum([0, 0, 1] * vec, axis=-1)
    threshold = np.cos(angle_threshold).item()
    obs = vec[dotp >= threshold]

    camera_params = ruststartracker.CameraParameters(
        np.array(
            (
                (4000, 0, 900),
                (0, 4000, 500),
                (0, 0, 1),
            ),
            dtype=np.float32,
        ),
        (1800, 1000),
        None,
    )

    pixel_in_frame = (camera_params.camera_matrix @ obs.T).T
    pixel_in_frame = pixel_in_frame[..., :2] / pixel_in_frame[..., 2:]

    size_x, size_y = camera_params.cam_resolution
    img = np.zeros((size_y, size_x), np.uint8)
    for x, y in pixel_in_frame.astype(int):
        image_patch = img[y - 1 : y + 2, x - 1 : x + 2]
        image_patch[:] = 50

    return img, vec, mag, pixel_in_frame, camera_params


def test_star_matcher_success(setup):
    img, vec, mag, _, camera_params = setup

    rot = scipy.spatial.transform.Rotation.from_rotvec([1, 1, 1])
    vec = rot.inv().apply(vec)

    st = ruststartracker.StarTracker(
        vec,
        mag,
        camera_params,
        inter_star_angle_tolerance=np.radians(0.1).item(),
        n_minimum_matches=6,
    )
    res = st.process_image(img)

    assert res is not None
    np.testing.assert_allclose(res.quat, rot.inv().as_quat(), rtol=0.001, atol=0.001)
    assert res.n_matches >= 4


def test_star_matcher_exhaust(setup):
    img, vec, mag, _, camera_params = setup
    st = ruststartracker.StarTracker(
        vec,
        mag,
        camera_params,
        inter_star_angle_tolerance=np.radians(0.001).item(),
        n_minimum_matches=500,
        timeout_secs=999.0,
    )
    with pytest.raises(ruststartracker.StarTrackerError, match="SearchExhausted"):
        st.process_image(img)


def test_star_matcher_timout(setup):
    img, vec, mag, _, camera_params = setup
    timeout = 0.0002
    st = ruststartracker.StarTracker(
        vec,
        mag,
        camera_params,
        inter_star_angle_tolerance=np.radians(0.1).item(),
        n_minimum_matches=500,
        timeout_secs=timeout,
    )
    t = time.monotonic()
    with pytest.raises(ruststartracker.StarTrackerError, match="Timeout"):
        st.process_image(img)
    passed_time = time.monotonic() - t
    assert passed_time > timeout


def test_star_matcher_not_enough_stars(setup):
    _, vec, mag, pixel_in_frame, camera_params = setup
    timeout = 0.2
    st = ruststartracker.StarTracker(
        vec,
        mag,
        camera_params,
        inter_star_angle_tolerance=np.radians(0.1).item(),
        n_minimum_matches=500,
        timeout_secs=timeout,
    )
    with pytest.raises(ruststartracker.StarTrackerError, match="NotEnoughStars"):
        st.process_image_coordiantes(pixel_in_frame[:2])


if __name__ == "__main__":
    pytest.main([__file__, "--capture=no"])
