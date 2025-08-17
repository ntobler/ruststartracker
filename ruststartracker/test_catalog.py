import datetime
import time

import astropy.time  # type: ignore[import]
import numpy as np
import pytest

import ruststartracker.catalog
import ruststartracker.libruststartracker


def test_time_to_epoch():
    np.testing.assert_allclose(
        ruststartracker.catalog.time_to_epoch(
            datetime.datetime.fromisoformat("2000-01-01T11:58:56")
        ),
        2000.0,
        rtol=1e-20,
        atol=1 / (365 * 86400),
    )


@pytest.mark.parametrize(
    "iso_date",
    [
        "2000-01-01T11:58:56",
        "2024-01-01T11:58:56",
    ],
)
def test_time_to_epoch_astropy(iso_date: str):
    ground_truth = float(astropy.time.Time(iso_date).jyear)  # type: ignore
    np.testing.assert_allclose(
        ruststartracker.catalog.time_to_epoch(datetime.datetime.fromisoformat(iso_date)),
        ground_truth,
        rtol=1e-20,
        atol=64 / (365 * 86400),
    )


def test_extract_observations():
    positions = ruststartracker.catalog.StarCatalog.from_gaia().normalized_positions()

    assert positions.ndim == 2
    assert positions.shape[1] == 3
    np.testing.assert_allclose(np.linalg.norm(positions, axis=-1), 1.0, rtol=1e-5)


def test_gaia_python_rust():
    t0 = time.monotonic()
    positions = ruststartracker.catalog.StarCatalog.from_gaia().normalized_positions(epoch=2025.0)
    print(f"Python catalog took {time.monotonic() - t0:.3f} seconds")
    t0 = time.monotonic()
    positions2 = ruststartracker.libruststartracker.StarCatalog.from_gaia(
        max_magnitude=6.0
    ).normalized_positions(epoch=2025.0, observer_position=None)
    print(f"Rust catalog took {time.monotonic() - t0:.3f} seconds")
    np.testing.assert_allclose(positions, positions2, rtol=1e-5, atol=1e-5)


def test_hipparcos_python_rust():
    t0 = time.monotonic()
    positions = ruststartracker.catalog.StarCatalog.from_hipparcos().normalized_positions(
        epoch=2025.0
    )
    print(f"Python catalog took {time.monotonic() - t0:.3f} seconds")
    t0 = time.monotonic()
    positions2 = ruststartracker.libruststartracker.StarCatalog.from_hipparcos(
        max_magnitude=6.0
    ).normalized_positions(epoch=2025.0, observer_position=None)
    print(f"Rust catalog took {time.monotonic() - t0:.3f} seconds")
    np.testing.assert_allclose(positions, positions2, rtol=1e-5, atol=1e-5)


def test_compare_hipparcos_gaia():
    def filt(cat: ruststartracker.catalog.StarCatalog, m: float = 0):
        positions = cat.normalized_positions(epoch=2025.0)
        mags = cat.magnitude
        mask = (positions[..., 2] > 0.9) * (mags > 5 - m) * (mags < 8 - m)
        return positions[mask]

    positions_hipparcos = filt(ruststartracker.catalog.StarCatalog.from_hipparcos(max_magnitude=8))
    positions_gaia = filt(ruststartracker.catalog.StarCatalog.from_gaia(max_magnitude=8), m=0.9)

    dists = np.linalg.norm(
        positions_hipparcos[np.newaxis, :, :] - positions_gaia[:, np.newaxis, :], axis=-1
    )
    matches = np.min(dists, axis=0) < 0.00005

    assert np.mean(matches > 0.9)

    if False:
        import matplotlib.pyplot as plt

        plt.plot(positions_gaia[..., 0], positions_gaia[..., 1], "+", label="gaia")
        plt.plot(positions_hipparcos[..., 0], positions_hipparcos[..., 1], "x", label="hipparcos")
        plt.legend()
        plt.show()


if __name__ == "__main__":
    pytest.main([__file__])
