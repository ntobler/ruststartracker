"""Create star catalog."""

import csv
import datetime
import math
import pathlib

import numpy as np
import numpy.typing as npt

AU: float = 149597870.693
"""Astronomical unit."""

INTERNAL_CATALOG_FILE = (
    pathlib.Path(__file__).parent.expanduser().absolute() / "gaia_data_j2016.csv"
)
"""Location of the internal star catalog file."""


def time_to_epoch(t: datetime.datetime) -> float:
    """Convert a date time object to the Epoch in Julian years."""
    # Get time stamp at the Julian year J2000
    delta_t_j2000 = 64  # Delta T at J2000 https://en.wikipedia.org/wiki/%CE%94T_(timekeeping)
    j2000 = datetime.datetime.fromisoformat("2000-01-01T12:00:00").timestamp() - delta_t_j2000
    # Calculate difference from given datetime to j2000 and convert to Julian years
    seconds_in_j_year = 365.25 * 86400
    return (t.timestamp() - j2000) / seconds_in_j_year + 2000.0


class StarCatalog:
    """Star catalog from Gaia data."""

    _data: npt.NDArray[np.float32]
    """Underlying data array."""
    ra: npt.NDArray[np.float32]
    """Right ascension at epoch in rads."""
    de: npt.NDArray[np.float32]
    """Declination at epoch in rads."""
    parallax: npt.NDArray[np.float32]
    """Parallax in rads."""
    proper_motion_ra: npt.NDArray[np.float32]
    """Proper motion of the right ascension in mad."""
    proper_motion_de: npt.NDArray[np.float32]
    """Proper motion of the declination in mad."""
    magnitude: npt.NDArray[np.float32]
    """Magnitude values."""
    epoch: float = 2016.0
    """Epoch of the catalog in years."""

    def __init__(
        self, filename: pathlib.Path | str | None = None, max_magnitude: float = 6.0
    ) -> None:
        """Read catalog from file.

        Args:
            filename: Star catalog filename.
            max_magnitude: Maximum magnitude to include.
        """
        # Use locally included star catalog is no file is given
        filename = INTERNAL_CATALOG_FILE if filename is None else pathlib.Path(filename)

        keep_columns = ("ra", "dec", "parallax", "pmra", "pmdec", "phot_g_mean_mag")

        with filename.open("r") as f:
            it = csv.reader(f, delimiter=",", strict=True)

            # Get columns
            columns = next(it)

            keep_columns_indices = tuple(columns.index(x) for x in keep_columns)
            min_length = len(keep_columns_indices)

            mag_column_index = columns.index("phot_g_mean_mag")

            rows = [
                [float(line[j]) for j in keep_columns_indices]
                for line in it
                if len(line) >= min_length and float(line[mag_column_index]) <= max_magnitude
            ]
        self._data = np.array(rows, dtype=np.float32)

        self.ra = self._data[:, keep_columns.index("ra")]
        self.de = self._data[:, keep_columns.index("dec")]
        self.parallax = self._data[:, keep_columns.index("parallax")]
        self.proper_motion_ra = self._data[:, keep_columns.index("pmra")]
        self.proper_motion_de = self._data[:, keep_columns.index("pmdec")]
        self.magnitude = self._data[:, keep_columns.index("phot_g_mean_mag")]

        deg2rad = math.pi / 180
        arcsec2deg = 1 / 3600
        mas2arcsec = 1 / 1000
        mas2rad = mas2arcsec * arcsec2deg * deg2rad

        # Convert data to rad
        self.ra *= deg2rad
        self.de *= deg2rad
        self.proper_motion_ra *= mas2rad
        self.proper_motion_de *= mas2rad
        self.parallax *= mas2rad

    def normalized_positions(
        self, *, epoch: float | None = None, observer_position: np.ndarray | None = None
    ) -> npt.NDArray[np.float32]:
        """Get star positions as normalized (x, y, z) vector.

        Args:
            epoch: The Julian epoch at which the positions should be calculated in Julian
                years. E.g. 2024.3. If None, the current local time us used.
            observer_position: The position of the observer in the Equatorial coordinate
                system relative to the sun. If None, position is not corrected.

        Returns:
            Normalized (x, y, z) vector of star positions, shape=[n, 3]
        """
        if epoch is None:
            epoch = time_to_epoch(datetime.datetime.now(tz=datetime.timezone.utc))

        delta_epoch = epoch - self.epoch

        # Precalculate sin and cos
        cos_ra: npt.NDArray[np.float32] = np.cos(self.ra)
        sin_ra: npt.NDArray[np.float32] = np.sin(self.ra)
        sin_de: npt.NDArray[np.float32] = np.sin(self.de)
        cos_de: npt.NDArray[np.float32] = np.cos(self.de)
        zeros = np.zeros_like(self.de)

        # Get star positions as normalized vector (x, y, z)
        vectors: npt.NDArray[np.float32] = np.stack(
            [cos_de * cos_ra, cos_de * sin_ra, sin_de], axis=-1
        )

        # Correct proper motion
        p_hat = np.stack([-sin_ra, cos_ra, zeros], axis=-1)
        q_hat = np.stack([-sin_de * cos_ra, -sin_de * sin_ra, cos_de], axis=-1)
        pm = delta_epoch * (
            self.proper_motion_ra[..., np.newaxis] * p_hat
            + self.proper_motion_de[..., np.newaxis] * q_hat
        )
        vectors += pm

        if observer_position is not None:
            plx = self.parallax[:, np.newaxis] * (observer_position / AU)[np.newaxis, :]
            vectors -= plx

        # normalize star positions
        vectors /= np.linalg.norm(vectors, axis=-1, keepdims=True)
        return vectors
