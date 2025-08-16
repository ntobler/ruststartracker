"""Star tracker implementation."""

import dataclasses

import cv2
import numpy as np
import numpy.typing as npt

import ruststartracker.libruststartracker


class StarTrackerError(Exception):
    """Exception related to the StarTracker."""


@dataclasses.dataclass
class CameraParameters:
    """Camera calibration parameters."""

    camera_matrix: npt.NDArray[np.float32]
    """3x3 intrinsic camera matrix."""
    cam_resolution: tuple[int, int]
    """Image resolution width, height"""
    dist_coefs: npt.NDArray[np.float32] | None
    """OpenCV distortion coefficients."""
    camera_mat_inv: npt.NDArray[np.float32] = dataclasses.field(init=False)
    """Inverted 3x3 intrinsic camera matrix."""

    def __post_init__(self) -> None:
        """Set inferred attributes."""
        self.camera_mat_inv = np.linalg.inv(self.camera_matrix)


@dataclasses.dataclass(frozen=True)
class StarTrackerResult:
    """Result container containing attitude and other information."""

    quat: npt.NDArray[np.float32]
    """Quaternion attitude (i, j, k, w)"""
    match_ids: npt.NDArray[np.uint32]
    """IDs of the matched stars."""
    n_matches: int
    """Number of matched stars"""
    duration_s: float
    """Processing duration in seconds."""
    mached_obs_x: npt.NDArray[np.float32]
    """Matched observation coordinates (x, y, z) , shape=[n_matches, 3]."""
    obs_indices: npt.NDArray[np.uint32]
    """Observation point indices of matches."""


def _image_coords_to_normed_vectors(
    camera_params: CameraParameters, image_coords: npt.NDArray[np.float32]
) -> npt.NDArray[np.float32]:
    """Convert image coordinates to normalized object coordinates.

    Args:
        camera_params: Camera calibration parameters
        image_coords: 2-dimensional image coordinates, shape=[n, 2]

    Returns:
        Normalized image coordinates, shape=[n, 3]
    """
    # Correct corner points with lens distortion if available
    if camera_params.dist_coefs is not None:
        image_coords_corrected = cv2.undistortPoints(
            src=image_coords[:, np.newaxis],
            cameraMatrix=camera_params.camera_matrix,
            distCoeffs=camera_params.dist_coefs,
            P=camera_params.camera_matrix,
        )[:, 0]
    else:
        image_coords_corrected = image_coords

    # Convert to 3-dimensional unit vectors (+z is the camera direction)
    homogeneous = np.concatenate(
        (image_coords_corrected, np.ones_like(image_coords_corrected[..., :1])), axis=-1
    )
    corner_coords_xyz: npt.NDArray[np.float32] = (camera_params.camera_mat_inv @ homogeneous.T).T
    corner_coords_xyz /= np.linalg.norm(corner_coords_xyz, axis=-1, keepdims=True)
    return corner_coords_xyz


class StarTracker:
    """Startracker attitude estimation."""

    def __init__(
        self,
        stars_xyz: npt.NDArray[np.float32],
        stars_mag: npt.NDArray[np.float32],
        camera_params: CameraParameters,
        *,
        max_lookup_magnitude: float | None = None,
        max_inter_star_angle: float | None = None,
        inter_star_angle_tolerance: float = 0.0008,
        n_minimum_matches: int = 10,
        timeout_secs: float = 1.0,
    ) -> None:
        """Set up star tracker environment for faster computation.

        Args:
            stars_xyz: Positions of catalog stars
            stars_mag: Magnitudes of catalog stars
            camera_params: Calibrated camera parameters
            max_lookup_magnitude: Maximum magnitude of stars used in the triangulation. Reducing
                this number means only bright stars are used for triangulation. This results in
                faster lookup performance.
            max_inter_star_angle: Maximum angle between stars that should be indexed.
                Calculating large inter star angles is expensive. If None, the angle is
                calculated from the camera field of view.
            inter_star_angle_tolerance: Tolerance for inter star angle matching in rad.
            n_minimum_matches: Minimum amount of required matches for a successful
                attitude estimation
            timeout_secs: Maximum allowed search time in seconds. A StarTrackerError is raised
                when the timeout elapses
        """
        self._camera_params = camera_params
        self._camera_mat_inv = np.linalg.inv(self._camera_params.camera_matrix)

        if max_inter_star_angle is None:
            # Get image corner points
            w, h = self._camera_params.cam_resolution
            corner_coords_pix = np.array(
                [(0, 0), (0, h - 1), (w - 1, 0), (w - 1, h - 1)], dtype=np.float32
            )

            # Convert to vectors
            corner_coords_xyz = _image_coords_to_normed_vectors(
                self._camera_params, corner_coords_pix
            )

            # Max angle is twice the angle from center to the farthest corner
            dot_products = (corner_coords_xyz * np.array([0, 0, 1], dtype=np.float32)).sum(axis=-1)
            max_inter_star_angle = float(np.arccos(dot_products.max())) * 2

        if max_lookup_magnitude is None:
            max_lookup_magnitude = 100.0  # A very faint star. Almost infinity

        self._star_matcher = ruststartracker.libruststartracker.StarMatcher(
            np.ascontiguousarray(stars_xyz, dtype=np.float32),
            np.ascontiguousarray(stars_mag, dtype=np.float32),
            float(max_lookup_magnitude),
            float(max_inter_star_angle),
            float(inter_star_angle_tolerance),
            int(n_minimum_matches),
            float(timeout_secs),
        )

    def get_centroids(
        self,
        img: npt.NDArray[np.uint8],
        darkframe: npt.NDArray[np.uint8] | None = None,
        n_candidates: int = 30,
        threshold: int | None = None,
        min_star_area: int = 4,
        max_star_area: int = 36,
    ) -> npt.NDArray[np.float32]:
        """Get star centroid image coordinates from a camera image."""
        if (
            (img.ndim != 2)
            or (img.dtype != np.uint8)
            or ((img.shape[1], img.shape[0]) != self._camera_params.cam_resolution)
        ):
            raise ValueError(
                f"img must be an uint8 array with shape {self._camera_params.cam_resolution}, "
                "matching the camera parameters"
            )

        if darkframe is not None:
            if (
                (darkframe.ndim != 2)
                or (darkframe.dtype != np.uint8)
                or ((darkframe.shape[1], darkframe.shape[0]) != self._camera_params.cam_resolution)
            ):
                raise ValueError(
                    "darkframe must be an uint8 array with shape "
                    f"{self._camera_params.cam_resolution}, matching the camera parameters"
                )
            # also clips the image
            img = cv2.subtract(img, darkframe)  # type: ignore[assignment]

        if min_star_area >= max_star_area:
            raise ValueError("min_star_area must be less than max_star_area")

        centroids, intensities = _extract_observations(
            img,
            threshold=threshold,
            min_star_area=min_star_area,
            max_star_area=max_star_area,
        )

        if threshold is None:
            threshold = ruststartracker.libruststartracker.get_threshold_from_histogram(
                img, fraction=0.99
            )

        centroids, intensities = ruststartracker.libruststartracker.extract_observations(
            img,
            threshold,
            min_star_area,
            max_star_area,
        )
        centroids = np.array(centroids, dtype=np.float32)
        intensities = np.array(intensities, dtype=np.float32)

        # At least 3 observations are required (one triangle)
        if len(centroids) < 3:
            raise StarTrackerError("Found too few star candidates (< 3) to continue.")

        # Sort candidates by their intensity. Brighter candidates are more likely to be stars
        bright_star_idx = intensities.argsort()[::-1][:n_candidates]
        # Limit number of candidates
        return centroids[bright_star_idx, :]

    def process_image(
        self,
        img: npt.NDArray[np.uint8],
        darkframe: npt.NDArray[np.uint8] | None = None,
        n_candidates: int = 30,
        threshold: int | None = None,
        min_star_area: int = 4,
        max_star_area: int = 36,
    ) -> StarTrackerResult:
        """Estimate attitude given a camera image."""
        centroids = self.get_centroids(
            img, darkframe, n_candidates, threshold, min_star_area, max_star_area
        )
        # Convert image star centroids to 3-dimensional unit vectors (+z is the camera direction)
        x_obs = _image_coords_to_normed_vectors(self._camera_params, centroids)
        return self.process_observation_vectors(x_obs)

    def process_image_coordiantes(self, image_xy: npt.NDArray[np.floating]) -> StarTrackerResult:
        """Estimate attitude given star image coordinates.

        Args:
            image_xy: Image positions of stars, shape=[n, 2].

        Returns:
            Result container with attitude.
        """
        x_obs = _image_coords_to_normed_vectors(self._camera_params, image_xy)
        return self.process_observation_vectors(x_obs)

    def process_observation_vectors(self, x_obs: npt.NDArray[np.floating]) -> StarTrackerResult:
        """Estimate attitude given star observations.

        Args:
            x_obs: Observation vectors, shape=[n, 3]. Must be normalized.

        Returns:
            Result container with attitude.
        """
        # Ensure input is c contiguous for rust backend.
        x_obs = np.ascontiguousarray(x_obs, dtype=np.float32)

        # Match stars in rust backend (numpy can't keep up with this)
        try:
            result = self._star_matcher.find(x_obs)
        except RuntimeError as e:
            # TODO provide better diagnostic e.g. minimum error etc.
            raise StarTrackerError(*e.args) from e

        quat, match_ids, obs_indices, n_matches, matched_obs, duration_s = result

        return StarTrackerResult(
            quat=quat,
            match_ids=match_ids,
            n_matches=n_matches,
            duration_s=duration_s,
            mached_obs_x=matched_obs,
            obs_indices=obs_indices,
        )


def _extract_observations(
    img: npt.NDArray[np.uint8],
    *,
    threshold: int | None = None,
    min_star_area: int = 3,
    max_star_area: int = 400,
) -> tuple[npt.NDArray[np.float32], npt.NDArray[np.float32]]:
    """Estimate sub pixel center and intensity from patches containing stars.

    Args:
        img: Grayscale image of stars
        threshold: Threshold level to find star patches. If None, it is chosen automatically.
        min_star_area: Minimum area of a patch, such that it is considered a star observation
        max_star_area: Maximum area of a patch, such that it is considered a star observation

    Returns:
        Sub-pixel star centers and an arbitrary metric for its brightness.
    """
    if threshold is None:
        # If threshold is not given, use 99th percentile from histogram
        hist = cv2.calcHist([img], [0], None, [256], [0, 256])[:, 0]
        cdf = hist.cumsum()
        cdf_normalized = cdf / cdf.max()
        threshold = np.argmax(cdf_normalized > 0.99).item()

    # convert the grayscale image to binary image
    _, img_binary = cv2.threshold(img, threshold, 1, cv2.THRESH_BINARY)

    # find connected components (patches of star candidates)
    _, _, stats, _ = cv2.connectedComponentsWithStats(img_binary, connectivity=4, ltype=cv2.CV_32S)
    # The first item is the background patch, which needs to be skipped
    stats = stats[1:]

    # Proceed with sub pixel center estimation
    # Only use regions that are within the area bounds
    areas = stats[:, cv2.CC_STAT_AREA]
    mask = areas >= min_star_area
    mask *= areas <= max_star_area
    stats = stats[mask]

    # initialize result arrays
    centers = np.empty((len(stats), 2), dtype=np.float32)
    intensities = np.empty((len(stats),), dtype=np.float32)

    for i, label in enumerate(stats):
        # Get region of interest (roi) bounding box
        top_pixel = label[cv2.CC_STAT_TOP]
        left_pixel = label[cv2.CC_STAT_LEFT]
        bottom_pixel = min(label[cv2.CC_STAT_HEIGHT] + top_pixel, img.shape[0] - 1)
        right_pixel = min(label[cv2.CC_STAT_WIDTH] + left_pixel, img.shape[1] - 1)

        # mask of intensities in ROI, in [row, col]
        roi = img[top_pixel : bottom_pixel + 1, left_pixel : right_pixel + 1]

        # if the image slice doesn't contain any pixels, continue
        if roi.size == 0:
            intensities[i] = 0
            continue

        intensity_sum = roi.sum(dtype=np.uint32).item()
        intensities[i] = intensity_sum

        x_bar = (np.arange(left_pixel, right_pixel + 1, dtype=np.uint32) * roi.sum(axis=0)).sum(
            dtype=np.float32
        ).item() / intensity_sum
        y_bar = (np.arange(top_pixel, bottom_pixel + 1, dtype=np.uint32) * roi.sum(axis=1)).sum(
            dtype=np.float32
        ).item() / intensity_sum
        centers[i] = (x_bar, y_bar)

    # Only return centers with intensity
    mask = intensities > 0
    return centers[mask], intensities[mask]
