import os

import numpy as np
import pytest
import scipy.spatial

import ruststartracker
import ruststartracker.star


@pytest.fixture()
def prepare() -> tuple[ruststartracker.StarTracker, np.ndarray]:
    os.environ["RUST_BACKTRACE"] = "1"

    camera_matrix = np.array(
        [
            [2.12694338e03, 0.00000000e00, 4.71566027e02],
            [0.00000000e00, 2.12548786e03, 3.10739295e02],
            [0.00000000e00, 0.00000000e00, 1.00000000e00],
        ]
    )

    cam_resolution = (960, 540)

    dist_coefs = np.array([-0.44120807, -0.15954202, 0.00767012, -0.00213292, -1.64788247])

    catalog = ruststartracker.StarCatalog()
    star_catalog_vecs = catalog.normalized_positions()

    camera_params = ruststartracker.CameraParameters(
        camera_matrix=camera_matrix,
        cam_resolution=cam_resolution,
        dist_coefs=dist_coefs,
    )

    st = ruststartracker.StarTracker(
        star_catalog_vecs,
        camera_params,
        inter_star_angle_tolerance=np.radians(0.1).item(),
        n_minimum_matches=5,
    )

    return st, star_catalog_vecs


def test_example(prepare: tuple[ruststartracker.StarTracker, np.ndarray]):
    os.environ["RUST_BACKTRACE"] = "1"

    st, _ = prepare

    obs = np.array(
        [
            [0.11975033, -0.02227603, 0.9925541],
            [0.03917335, 0.04533212, 0.99820361],
            [0.05137746, -0.01717139, 0.99853167],
            [-0.14742009, 0.00734109, 0.98904673],
            [0.03396359, 0.05851033, 0.99770888],
            [-0.10286126, 0.04479652, 0.99368649],
            [-0.050927, -0.06002669, 0.9968968],
            [0.02815389, 0.02852981, 0.99919638],
            [-0.1390861, 0.05684676, 0.98864731],
            [-0.14276463, -0.13402131, 0.98064089],
        ]
    )

    result = st.process_observation_vectors(obs)
    np.testing.assert_allclose(
        result.quat, [0.1722, -0.4309, 0.8818, 0.08396], rtol=1e-3, atol=1e-3
    )


def test_star_matcher(prepare: tuple[ruststartracker.StarTracker, np.ndarray]):
    os.environ["RUST_BACKTRACE"] = "1"

    st, vec = prepare

    angle_threshold = np.radians(10)
    dotp = np.sum([0, 0, 1] * vec, axis=-1)
    threshold = np.cos(angle_threshold).item()
    obs = vec[dotp >= threshold]

    rot = scipy.spatial.transform.Rotation.from_rotvec([1, 1, 1])

    obs = rot.apply(obs)

    res = st.process_observation_vectors(obs)

    assert res is not None
    np.testing.assert_allclose(res.quat, rot.inv().as_quat(), rtol=0.001, atol=0.001)
    assert res.n_matches >= 4


if __name__ == "__main__":
    pytest.main([__file__])
