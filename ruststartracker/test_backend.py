import os

import numpy as np
import pytest
import scipy.spatial

from ruststartracker import libruststartracker


def test_triangle_finder():
    ab = np.array([[234, 5643], [1, 2], [2, 4], [3, 9], [2, 6]])
    ac = np.array([[345, 2343], [8, 2], [3, 4], [1, 7], [0, 5], [3, 1]])
    bc = np.array([[435, 4355], [1, 0], [4, 8], [8, 1], [1, 9]])

    f = libruststartracker.TriangleFinder(ab, ac, bc)
    assert f.get() == [1, 2, 8]
    assert list(libruststartracker.IterTriangleFinder(ab, ac, bc)) == [
        [1, 2, 8],
        [2, 4, 8],
        [3, 9, 1],
    ]


def test_unit_vector_lookup():
    rng = np.random.default_rng(42)
    n_vecs = 2617

    vec = rng.normal(size=[n_vecs, 3]).astype(np.float32)
    vec /= np.linalg.norm(vec, axis=-1, keepdims=True)

    uvl = libruststartracker.UnitVectorLookup(vec)

    keys = rng.normal(size=[10, 3]).astype(np.float32)
    keys /= np.linalg.norm(keys, axis=-1, keepdims=True)

    for key in keys:
        np.testing.assert_array_equal(
            uvl.lookup_nearest(key),
            np.linalg.norm(vec - key, axis=-1).argmin().item(),
        )

    angle_threshold = np.radians(15)

    threshold = np.cos(angle_threshold).item()
    results = []
    angles = []
    for a in range(len(vec)):
        dotp = np.sum(vec[a] * vec[a + 1 :], axis=-1)
        b = np.nonzero(dotp >= threshold)[0]
        results.append(np.array([np.full(len(b), a), (a + 1) + b]))
        angles.append(np.arccos(dotp[b]))
    angles = np.concatenate(angles, axis=0)
    args = np.argsort(angles)
    close_indices_gt = np.concatenate(results, axis=-1).T[args]
    angles_gt = angles[args]

    close_indices, angles, poly = uvl.get_inter_star_index(
        np.array(vec[:, :3], dtype=np.float32),
        np.ones(len(vec), dtype=np.float32),
        angle_threshold,
        10,
    )
    close_indices = np.array(close_indices)
    angles = np.array(angles)
    poly = np.array(poly)

    # There are some cases where the float32 accuracy is insufficient to tell
    # angles apart. Consequently the order may be alightly different. However,
    # we're able to test if the not-matching indices align with items that have
    # at minimum one other angle of the exact same value
    i = (close_indices != close_indices_gt).any(axis=-1)
    assert (np.unique(angles[i], return_counts=True)[-1] >= 2).all()

    np.testing.assert_allclose(angles, angles_gt)

    actual_index = np.arange(angles_gt.size)

    poly_gt = np.polyfit(angles_gt, actual_index, 2)
    lookup_index = np.polyval(poly_gt, angles_gt)

    poly_gt[-1] -= (lookup_index - actual_index).max()
    lookup_index = np.polyval(poly_gt, angles_gt)

    if False:
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(2)
        axs[0].plot(angles, lookup_index, "--")
        axs[0].plot(angles, actual_index)
        axs[1].plot(lookup_index - actual_index)
        fig.savefig("angles.png")
        plt.close(fig)

    np.testing.assert_allclose(poly, poly_gt[::-1], rtol=0.001)


def test_star_matcher():
    rng = np.random.default_rng(42)

    os.environ["RUST_BACKTRACE"] = "1"

    n_cat_stars = 2617

    vec = rng.normal(size=[n_cat_stars, 3]).astype(np.float32)
    vec /= np.linalg.norm(vec, axis=-1, keepdims=True)

    magnitudes = rng.uniform(0, 10, size=vec.shape[:1]).astype(np.float32)

    key = rng.normal(size=[3]).astype(np.float32)
    key /= np.linalg.norm(key, axis=-1, keepdims=True)

    angle_threshold = np.radians(7)
    dotp = np.sum(key * vec, axis=-1)
    threshold = np.cos(angle_threshold).item()
    b = np.nonzero(dotp >= threshold)[0]
    obs_index = rng.permutation(b)
    obs = vec[obs_index]

    rot = scipy.spatial.transform.Rotation.from_rotvec([1, 1, 1])

    obs_rotated = rot.apply(obs).astype(np.float32)

    index = libruststartracker.StarMatcher(
        vec,
        magnitudes,
        10,
        np.radians(10).item(),
        np.radians(0.1).item(),
        4,
        999.0,
    )

    res = index.find(obs_rotated)

    assert res is not None

    quat, match_ids, obs_indices, n_matches, matched_obs, time_s = res
    np.testing.assert_allclose(quat, rot.inv().as_quat(), rtol=1e-6)
    assert n_matches >= 4
    assert len(obs_index) == len(match_ids)


if __name__ == "__main__":
    pytest.main([__file__])
