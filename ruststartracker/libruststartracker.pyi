from collections.abc import Iterator

import numpy as np
import numpy.typing as npt
from typing_extensions import Self

class StarMatcher:
    def __init__(
        self,
        stars_xyz: npt.NDArray[np.float32],
        stars_mag: npt.NDArray[np.float32],
        max_inter_star_angle: float,
        max_lookup_magnitude: float,
        inter_star_angle_tolerance: float,
        n_minimum_matches: int,
        timeout_secs: float,
    ) -> None: ...
    def find(
        self, obs_xyz: npt.NDArray[np.float32]
    ) -> tuple[
        npt.NDArray[np.float32],
        npt.NDArray[np.uint32],
        npt.NDArray[np.uint32],
        int,
        npt.NDArray[np.float32],
        float,
    ]: ...

class TriangleFinder:
    def __init__(
        self,
        ab: npt.NDArray[np.float32],
        ac: npt.NDArray[np.float32],
        bc: npt.NDArray[np.float32],
    ) -> None: ...
    def get(self) -> list[int]: ...

class IterTriangleFinder:
    def __init__(
        self,
        ab: npt.NDArray[np.float32],
        ac: npt.NDArray[np.float32],
        bc: npt.NDArray[np.float32],
    ) -> None: ...
    def __iter__(self) -> Iterator[list[int]]: ...

class UnitVectorLookup:
    def __init__(self, vec: npt.NDArray[np.float32]) -> None: ...
    def lookup_nearest(self, key: npt.NDArray[np.float32]) -> int: ...
    def get_inter_star_index(
        self,
        stars: npt.NDArray[np.float32],
        magnitudes: npt.NDArray[np.float32],
        max_angle_rad: float,
        max_magnitude: float,
    ) -> tuple[list[list[int]], list[float], list[float]]: ...
    def look_up_close_angles(
        self,
        vectors: npt.NDArray[np.float32],
        magnitudes: npt.NDArray[np.float32],
        max_angle_rad: float,
        max_magnitude: float,
    ) -> list[tuple[list[float], float]]: ...
    def look_up_close_angles_naive(
        self,
        vectors: npt.NDArray[np.float32],
        magnitudes: npt.NDArray[np.float32],
        max_angle_rad: float,
        max_magnitude: float,
    ) -> list[tuple[list[float], float]]: ...

def get_threshold_from_histogram(
    img: npt.NDArray[np.uint8],
    *,
    fraction: float,
) -> int: ...
def extract_observations(
    img: npt.NDArray[np.uint8],
    threshold: int,
    min_star_area: int,
    max_star_area: int,
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.float32]]: ...

class StarCatalog:
    @classmethod
    def from_gaia(cls, *, max_magnitude: float | None) -> Self: ...
    def normalized_positions(
        self, *, epoch: float | None, observer_position: np.ndarray | None
    ) -> npt.NDArray[np.float32]: ...
