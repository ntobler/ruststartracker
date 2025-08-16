use numpy::{self, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;
use pyo3::types::PyAny;
#[cfg(feature = "gaia")]
use pyo3::types::PyType;
use pyo3::{
    exceptions::PyRuntimeError, pyclass, pymethods, pymodule, types::PyModule, Bound, PyRef,
    PyRefMut, PyResult,
};
use std::path::PathBuf;
use std::{time::Instant, usize};

mod ordered_combinations;
pub mod star;
pub mod starcat;
#[cfg(feature = "improc")]
pub mod starextraction;
mod tree;
mod trianglefinder;

#[pymodule]
fn libruststartracker(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<TriangleFinder>()?;
    m.add_class::<IterTriangleFinder>()?;
    m.add_class::<StarMatcher>()?;
    m.add_class::<UnitVectorLookup>()?;
    m.add_class::<StarCatalog>()?;
    #[cfg(feature = "improc")]
    m.add_function(wrap_pyfunction!(get_threshold_from_histogram, m)?)?;
    #[cfg(feature = "improc")]
    m.add_function(wrap_pyfunction!(extract_observations, m)?)?;
    Ok(())
}

#[pyclass]
struct TriangleFinder {
    inner: trianglefinder::TriangleFinder,
}

#[pymethods]
impl TriangleFinder {
    #[new]
    fn new(
        connections_ab: Vec<[u32; 2]>,
        connections_ac: Vec<[u32; 2]>,
        connections_bc: Vec<[u32; 2]>,
    ) -> Self {
        TriangleFinder {
            inner: trianglefinder::TriangleFinder::new(
                connections_ab,
                connections_ac,
                connections_bc,
            ),
        }
    }

    pub fn get(&self) -> PyResult<Option<[u32; 3]>> {
        Ok(self.inner.get())
    }

    pub fn get_all(&self) -> PyResult<Vec<[u32; 3]>> {
        Ok(self.inner.get_all())
    }
}

#[pyclass]
struct IterTriangleFinder {
    iter: trianglefinder::IterTriangleFinder,
}

#[pymethods]
impl IterTriangleFinder {
    #[new]
    fn new(
        connections_ab: Vec<[u32; 2]>,
        connections_ac: Vec<[u32; 2]>,
        connections_bc: Vec<[u32; 2]>,
    ) -> Self {
        IterTriangleFinder {
            iter: trianglefinder::IterTriangleFinder::new(trianglefinder::TriangleFinder::new(
                connections_ab,
                connections_ac,
                connections_bc,
            )),
        }
    }

    fn __iter__(slf: PyRef<'_, Self>) -> PyRef<'_, Self> {
        slf
    }
    fn __next__(mut slf: PyRefMut<'_, Self>) -> Option<[u32; 3]> {
        slf.iter.next()
    }
}

#[pyclass]
struct StarMatcher {
    inner: star::StarMatcher,
}

#[pymethods]
impl StarMatcher {
    #[new]
    fn new<'py>(
        stars_xyz: numpy::PyReadonlyArray2<'py, f32>,
        stars_mag: numpy::PyReadonlyArray1<'py, f32>,
        max_lookup_magnitude: f32,
        max_inter_star_angle: f32,
        inter_star_angle_tolerance: f32,
        n_minimum_matches: usize,
        timeout_secs: f32,
    ) -> PyResult<Self> {
        let stars_slice: &[[f32; 3]] = numpy_to_slice_2d(&stars_xyz)?;
        let mags_slice: &[f32] = numpy_to_slice_1d(&stars_mag)?;
        match star::StarMatcher::new(
            stars_slice.to_vec(),
            mags_slice,
            max_lookup_magnitude,
            max_inter_star_angle,
            inter_star_angle_tolerance,
            n_minimum_matches,
            timeout_secs,
        ) {
            Ok(inner) => Ok(StarMatcher { inner }),
            Err(e) => Err(PyRuntimeError::new_err(e)),
        }
    }

    pub fn find<'py>(
        &self,
        py: Python<'py>,
        obs_xyz: numpy::PyReadonlyArray2<'py, f32>,
    ) -> PyResult<(
        Bound<'py, numpy::PyArray1<f32>>,
        Bound<'py, numpy::PyArray1<u32>>,
        Bound<'py, numpy::PyArray1<u32>>,
        u32,
        Bound<'py, numpy::PyArray2<f32>>,
        f32,
    )> {
        let obs_xyz_slice: &[[f32; 3]] = numpy_to_slice_2d(&obs_xyz)?;
        let now = Instant::now();
        let res = self
            .inner
            .find(obs_xyz_slice)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        let duration_s = now.elapsed().as_secs_f32();
        Ok((
            numpy_from_slice_1d(py, &res.quat),
            numpy_from_slice_1d(py, &res.match_ids),
            numpy_from_slice_1d(py, &res.obs_indices),
            res.n_matches,
            numpy_from_slice_2d(py, &res.obs_matched),
            duration_s,
        ))
    }
}

#[pyclass]
struct UnitVectorLookup {
    inner: tree::UnitVectorLookup,
}

#[pymethods]
impl UnitVectorLookup {
    #[new]
    fn new(vectors: Vec<[f32; 3]>) -> Self {
        UnitVectorLookup {
            inner: tree::UnitVectorLookup::new(&vectors),
        }
    }

    pub fn lookup_nearest(&self, vector: [f32; 3]) -> PyResult<usize> {
        Ok(self.inner.lookup_nearest(&vector))
    }

    pub fn get_inter_star_index<'py>(
        &self,
        stars: numpy::PyReadonlyArray2<'py, f32>,
        magnitudes: numpy::PyReadonlyArray1<'py, f32>,
        max_angle_rad: f32,
        max_magnitude: f32,
    ) -> PyResult<(Vec<[u32; 2]>, Vec<f32>, [f32; 3])> {
        let stars_slice: &[[f32; 3]] = numpy_to_slice_2d(&stars)?;
        let magnitudes_slice: &[f32] = numpy_to_slice_1d(&magnitudes)?;
        let now = Instant::now();
        let res = star::InterStarIndex::new(
            &self.inner,
            stars_slice,
            magnitudes_slice,
            max_angle_rad,
            max_magnitude,
        )
        .map_err(|e| {
            PyRuntimeError::new_err(format!("Could not calculate inter star angle: {}", e))
        })?;

        println!("Time passed: {:?}", now.elapsed());
        Ok((res.pairs, res.angles, res.polynomial))
    }

    pub fn look_up_close_angles(
        &self,
        vectors: Vec<[f32; 3]>,
        magnitudes: Vec<f32>,
        max_angle_rad: f32,
        max_magnitude: f32,
    ) -> PyResult<Vec<([u32; 2], f32)>> {
        let now = Instant::now();
        let res =
            self.inner
                .look_up_close_angles(&vectors, &magnitudes, max_angle_rad, max_magnitude);
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }

    pub fn look_up_close_angles_naive(
        &self,
        vectors: Vec<[f32; 3]>,
        magnitudes: Vec<f32>,
        max_angle_rad: f32,
        max_magnitude: f32,
    ) -> PyResult<Vec<([u32; 2], f32)>> {
        let now = Instant::now();
        let res =
            star::look_up_close_angles_naive(&vectors, &magnitudes, max_angle_rad, max_magnitude);
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }
}

#[pyclass]
struct StarCatalog {
    inner: starcat::StarCatalog,
}

#[pymethods]
impl StarCatalog {
    #[new]
    fn new(filename: Bound<'_, PyAny>, epoch: f64, max_magnitude: Option<f64>) -> PyResult<Self> {
        let path: PathBuf = filename.extract()?;
        Ok(StarCatalog {
            inner: starcat::StarCatalog::new_from_file(path, epoch, max_magnitude)
                .map_err(PyRuntimeError::new_err)?,
        })
    }

    #[cfg(feature = "gaia")]
    #[classmethod]
    fn from_gaia(_cls: &Bound<'_, PyType>, max_magnitude: Option<f64>) -> PyResult<Self> {
        Ok(StarCatalog {
            inner: starcat::StarCatalog::new_from_gaia(max_magnitude)
                .map_err(PyRuntimeError::new_err)?,
        })
    }

    pub fn normalized_positions(
        &self,
        epoch: Option<f64>,
        observer_position: Option<[f64; 3]>,
    ) -> Vec<[f64; 3]> {
        self.inner.normalized_positions(epoch, observer_position)
    }
}

#[cfg(feature = "improc")]
#[pyfunction]
pub fn get_threshold_from_histogram<'py>(
    img: numpy::PyReadonlyArray2<'py, u8>,
    fraction: f64,
) -> PyResult<u8> {
    if !img.is_c_contiguous() {
        return Err(PyRuntimeError::new_err("Image must be a c_contiguous"));
    }
    Ok(starextraction::get_threshold_from_histogram(
        img.as_slice()?,
        fraction,
    ))
}

#[cfg(feature = "improc")]
#[pyfunction]
pub fn extract_observations<'py>(
    py: Python<'py>,
    img: numpy::PyReadonlyArray2<'py, u8>,
    threshold_value: u8,
    min_area: usize,
    max_area: usize,
) -> PyResult<(
    Bound<'py, numpy::PyArray2<f64>>,
    Bound<'py, numpy::PyArray1<f64>>,
)> {
    if !img.is_c_contiguous() || img.ndim() != 2 {
        return Err(PyRuntimeError::new_err(
            "Image must be a c_contiguous 2D array",
        ));
    }
    let (centroids, intensities) = starextraction::extract_observations(
        img.as_slice()?,
        (img.shape()[1], img.shape()[0]),
        threshold_value,
        min_area,
        max_area,
    )
    .map_err(|e| PyRuntimeError::new_err(e))?;

    let centroids_np = numpy_from_slice_2d(py, &centroids);
    let intensities_np = numpy_from_slice_1d(py, &intensities);

    Ok((centroids_np, intensities_np))
}

fn numpy_from_slice_2d<'py, const N: usize, T>(
    py: Python<'py>,
    slice: &[[T; N]],
) -> Bound<'py, numpy::PyArray2<T>>
where
    T: numpy::Element,
{
    let len = slice.len(); // n
    let ptr = slice.as_ptr() as *const T;
    let total_len = len * N;
    let contiguous_slice = unsafe { std::slice::from_raw_parts(ptr, total_len) };
    let flat_np = contiguous_slice.to_pyarray_bound(py);
    flat_np.reshape((len, N)).unwrap() // Save to unwrap as we know the shape is correct
}

fn numpy_from_slice_1d<'py, T>(py: Python<'py>, slice: &[T]) -> Bound<'py, numpy::PyArray1<T>>
where
    T: numpy::Element,
{
    slice.to_pyarray_bound(py)
}

fn numpy_to_slice_2d<'py, T: numpy::Element + Copy, const L: usize>(
    array: &'py numpy::PyReadonlyArray2<'py, T>,
) -> PyResult<&'py [[T; L]]> {
    if !array.is_c_contiguous() || array.ndim() != 2 || array.shape()[1] != L {
        return Err(PyRuntimeError::new_err(format!(
            "vectors must be a c_contiguous array with shape=[n, {}]",
            L
        )));
    }
    let slice = array.as_slice()?;
    Ok(unsafe { std::slice::from_raw_parts(slice.as_ptr() as *const [T; L], slice.len() / L) })
}

fn numpy_to_slice_1d<'py, T: numpy::Element + Copy>(
    array: &'py numpy::PyReadonlyArray1<'py, T>,
) -> PyResult<&'py [T]> {
    if !array.is_c_contiguous() || array.ndim() != 1 {
        return Err(PyRuntimeError::new_err(
            "vectors must be a c_contiguous array with shape=[n]",
        ));
    }
    Ok(array.as_slice()?)
}
