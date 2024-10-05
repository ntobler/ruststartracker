import datetime

import astropy.time
import numpy as np
import pytest

import ruststartracker.catalog


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
    positions = ruststartracker.catalog.StarCatalog().normalized_positions()

    assert positions.ndim == 2
    assert positions.shape[1] == 3
    np.testing.assert_allclose(np.linalg.norm(positions, axis=-1), 1.0, rtol=1e-5)


if __name__ == "__main__":
    pytest.main([__file__])
