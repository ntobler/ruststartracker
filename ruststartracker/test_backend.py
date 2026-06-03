import os

import numpy as np
import pytest
import scipy.spatial

from ruststartracker import libruststartracker


def test_lookup_nearest():
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


def test_close_angle_lookup():
    rng = np.random.default_rng(42)
    n_vecs = 2617

    vec = rng.normal(size=[n_vecs, 3]).astype(np.float32)
    vec /= np.linalg.norm(vec, axis=-1, keepdims=True)

    uvl = libruststartracker.UnitVectorLookup(vec)

    angle_threshold = np.radians(15)

    threshold = np.cos(angle_threshold).item()
    pairs_gt = []
    angles_gt = []
    for a in range(len(vec)):
        dotp = np.sum(vec[a] * vec[a + 1 :], axis=-1)
        b = np.nonzero(dotp >= threshold)[0]
        pairs_gt.append(np.array([np.full(len(b), a), (a + 1) + b]))
        angles_gt.append(np.arccos(dotp[b]))
    angles_gt = np.concatenate(angles_gt, axis=0)
    pairs_gt = np.concatenate(pairs_gt, axis=-1).T
    args = np.argsort(angles_gt)
    pairs_gt = pairs_gt[args]
    angles_gt = angles_gt[args]

    res = uvl.look_up_close_angles(
        np.array(vec, dtype=np.float32),
        np.ones(len(vec), dtype=np.float32),
        np.cos(angle_threshold).item(),
        10,
    )
    pairs = np.array([r[0] for r in res])
    angles_proxy = np.array([r[1] for r in res])
    angles = np.arccos(angles_proxy)
    args = np.argsort(angles)
    pairs = pairs[args]
    angles = angles[args]

    # Check if the same angles are returned (order has already been normalized by the sorting)
    np.testing.assert_allclose(angles, angles_gt)

    # There are some cases where the float32 accuracy is insufficient to tell
    # angles apart. Consequently the order may be slightly different. However,
    # we're able to test if the not-matching indices align with items that have
    # at minimum one other angle of the exact same value
    i = (pairs != pairs_gt).any(axis=-1)
    assert np.mean(i) < 0.1, "More than 10 percent of pairs do not match"
    assert (np.unique(angles[i], return_counts=True)[-1] >= 2).all()


def test_get_inter_star_index():
    os.environ["RUST_BACKTRACE"] = "1"

    rng = np.random.default_rng(42)
    n_vecs = 2617

    vec = rng.normal(size=[n_vecs, 3]).astype(np.float32)
    vec /= np.linalg.norm(vec, axis=-1, keepdims=True)

    uvl = libruststartracker.UnitVectorLookup(vec)

    angle_threshold = np.radians(15)

    threshold = np.cos(angle_threshold).item()
    pairs_gt = []
    angles_gt = []
    for a in range(len(vec)):
        dotp = np.sum(vec[a] * vec[a + 1 :], axis=-1)
        b = np.nonzero(dotp >= threshold)[0]
        pairs_gt.append(np.array([np.full(len(b), a), (a + 1) + b]))
        angles_gt.append(np.arccos(dotp[b]))
    angles_gt = np.concatenate(angles_gt, axis=0)
    pairs_gt = np.concatenate(pairs_gt, axis=-1).T
    args = np.argsort(angles_gt)
    pairs_gt = pairs_gt[args]
    angles_gt = angles_gt[args]

    angle_proxy_gt = np.cos(angles_gt)

    lookup_center = np.radians(7).item()
    lookup_tolerance = np.radians(0.1).item()

    pairs, angles_proxy, poly, looked_up_pairs = uvl.get_inter_star_index(
        np.array(vec, dtype=np.float32),
        np.ones(len(vec), dtype=np.float32),
        angle_threshold,
        10,
        lookup_center,
        lookup_tolerance,
    )
    pairs = np.array(pairs)
    angles_proxy = np.array(angles_proxy)
    poly = np.array(poly)[::-1]  # Reverse polynomial to match numpy's polyval order
    looked_up_pairs = np.array(looked_up_pairs)

    # angles are monotonically increasing, so the angle proxy should decreasing
    assert np.diff(angles_proxy).max() <= 0, "Angle proxy is not sorted"
    np.testing.assert_allclose(angles_proxy, angle_proxy_gt, rtol=1e-4, atol=1e-12)

    # There are some cases where the float32 accuracy is insufficient to tell
    # angles apart. Consequently the order may be slightly different. However,
    # we're able to test if the not-matching indices align with items that have
    # at minimum one other angle of the exact same value
    i = (pairs != pairs_gt).any(axis=-1)
    assert np.mean(i) < 0.1, "More than 10 percent of pairs do not match"
    assert (np.unique(angles_proxy[i], return_counts=True)[-1] >= 2).all()

    scale = len(pairs_gt) + 1

    lookup_index = np.polyval(poly, (1.0 - angles_proxy) * 10.0) * scale

    transformed_angle_proxy_gt = (1.0 - angle_proxy_gt) * 10.0

    indices = np.arange(angle_proxy_gt.size)
    scaled_indices = indices * (1.0 / scale)
    poly_gt = np.polyfit(transformed_angle_proxy_gt, scaled_indices, 2)
    lookup_index_gt = np.polyval(poly_gt, transformed_angle_proxy_gt) * scale
    max_val = (lookup_index_gt - indices).max()
    min_val = (lookup_index_gt - indices).min()
    print(max_val, min_val)
    poly_gt[-1] -= max_val / scale
    lookup_index_gt = np.polyval(poly_gt, transformed_angle_proxy_gt) * scale

    if False:
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots()
        axs.plot(angle_proxy_gt, label="angle proxy gt")
        axs.plot(angles_proxy, label="angle proxy")
        axs.legend()

        fig, axs = plt.subplots(2, sharex=True)
        axs[0].plot(angle_proxy_gt, lookup_index_gt, "--", label="lookup index gt")
        axs[0].plot(angle_proxy_gt, indices, label="actual index gt")
        axs[0].plot(angles_proxy, lookup_index, "--", label="lookup index")
        axs[0].plot(angles_proxy, indices, label="actual index")
        axs[0].legend()

        axs[1].plot(angle_proxy_gt, lookup_index_gt - indices, label="lookup error gt")
        axs[1].plot(angles_proxy, lookup_index - indices, label="lookup error")
        axs[1].legend()

        fig, axs = plt.subplots()
        axs.plot(poly_gt, "x", label="polynomial fit gt")
        axs.plot(poly, "x", label="polynomial fit")
        axs.legend()
        plt.show()

    np.testing.assert_allclose(lookup_index_gt - indices, lookup_index - indices, atol=0.1)

    # Fail if the polynomial fit is not accurate enough to represent the angles closely
    assert (
        max_val - min_val < 400
    ), "Polynomial fit is not accurate enough to preserve order of angles"

    # Check if the polynomial fit is correct
    np.testing.assert_allclose(poly, poly_gt, rtol=1e-3, atol=1e-5)

    # Check if the looked up pairs are correct

    mask = (angles_gt > lookup_center - lookup_tolerance) * (
        angles_gt < lookup_center + lookup_tolerance
    )

    matching_pairs_gt = pairs_gt[mask]

    # Due to f32 precision limitations and cosines close to 1.0,
    # we may miss some pairs that are very close to the lookup center.
    min_len = min(len(matching_pairs_gt), len(looked_up_pairs))
    assert (
        min_len / len(matching_pairs_gt) > 0.95
    ), "looked up less than 95 percent of the correct pairs"

    # There are some cases where the float32 accuracy is insufficient to tell
    # angles apart. Consequently the order may be slightly different. However,
    # we're able to test if the not-matching indices align with items that have
    # at minimum one other angle of the exact same value
    i = (matching_pairs_gt[:min_len] != looked_up_pairs[:min_len]).any(axis=-1)
    assert np.mean(i) < 0.1, "More than 10 percent of pairs do not match"


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
