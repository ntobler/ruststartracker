from typing import Iterator

import numpy as np
import numpy.typing as npt

def main_scratch() -> None: ...

class StarMatcher:
    def __init__(
        self,
        stars_xyz: npt.NDArray[np.float32],
        max_inter_star_angle: float,
        inter_star_angle_tolerance: float,
        n_minimum_matches: int,
        timeout_secs: float,
    ) -> None: ...
    def find(
        self, obs_xyz: npt.NDArray[np.float32]
    ) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.uint32], int, list[list[float]], float]: ...

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
    def get_inter_star_index_numpy(
        self, vec: npt.NDArray[np.float32], angle_threshold: float
    ) -> tuple[list[list[int]], list[float], list[float]]: ...
    def get_inter_star_index(
        self, vec: npt.NDArray[np.float32], angle_threshold: float
    ) -> tuple[list[list[int]], list[float], list[float]]: ...
    def look_up_close_angles(
        self, vectors: npt.NDArray[np.float32], max_angle_rad: float
    ) -> list[tuple[list[float], float]]: ...
    def look_up_close_angles_naive(
        self, vectors: npt.NDArray[np.float32], max_angle_rad: float
    ) -> list[tuple[list[float], float]]: ...
