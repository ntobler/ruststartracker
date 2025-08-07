use numpy::{self, PyArrayMethods, PyUntypedArrayMethods, ToPyArray};
use pyo3::prelude::*;
use pyo3::types::PyAny;
use pyo3::{
    exceptions::PyRuntimeError, pyclass, pymethods, pymodule, types::PyModule, types::PyType,
    Bound, PyRef, PyRefMut, PyResult,
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
mod util;

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
    fn new(
        stars_xyz: Vec<[f32; 3]>,
        max_inter_star_angle: f32,
        inter_star_angle_tolerance: f32,
        n_minimum_matches: usize,
        timeout_secs: f32,
    ) -> PyResult<Self> {
        match star::StarMatcher::new(
            stars_xyz,
            max_inter_star_angle,
            inter_star_angle_tolerance,
            n_minimum_matches,
            timeout_secs,
        ) {
            Ok(inner) => Ok(StarMatcher { inner }),
            Err(e) => Err(PyRuntimeError::new_err(e)),
        }
    }

    pub fn find(
        &self,
        obs_xyz: Vec<[f32; 3]>,
    ) -> PyResult<([f32; 4], Vec<u32>, Vec<u32>, u32, Vec<[f32; 3]>, f32)> {
        let now = Instant::now();
        let res = self.inner.find(&obs_xyz);
        let duration_s = now.elapsed().as_secs_f32();
        match res {
            Err(x) => Err(PyRuntimeError::new_err(x.to_string())),
            Ok(x) => Ok((
                x.quat,
                x.match_ids,
                x.obs_indices,
                x.n_matches,
                x.obs_matched,
                duration_s,
            )),
        }
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

    pub fn get_inter_star_index(
        &self,
        vectors: Vec<[f32; 3]>,
        max_angle_rad: f32,
    ) -> PyResult<(Vec<[u32; 2]>, Vec<f32>, [f32; 3])> {
        let now = Instant::now();
        let res = match star::get_inter_star_index(&self.inner, &vectors, max_angle_rad) {
            Ok(res) => res,
            Err(s) => {
                return Err(PyRuntimeError::new_err(format!(
                    "Could not calculate inter star angle: {}",
                    s
                )))
            }
        };
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }

    pub fn get_inter_star_index_numpy<'py>(
        &self,
        vectors: numpy::PyReadonlyArray2<'py, f32>,
        max_angle_rad: f32,
    ) -> PyResult<(Vec<[u32; 2]>, Vec<f32>, [f32; 3])> {
        let now = Instant::now();
        let vectors_inner = numpy_to_vec_3_32f(&vectors).unwrap();
        let res = match star::get_inter_star_index(&self.inner, vectors_inner, max_angle_rad) {
            Ok(res) => res,
            Err(s) => {
                return Err(PyRuntimeError::new_err(format!(
                    "Could not calculate inter star angle: {}",
                    s
                )))
            }
        };
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }

    pub fn look_up_close_angles(
        &self,
        vectors: Vec<[f32; 3]>,
        max_angle_rad: f32,
    ) -> PyResult<Vec<([u32; 2], f32)>> {
        let now = Instant::now();
        let res = self.inner.look_up_close_angles(&vectors, max_angle_rad);
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }

    pub fn look_up_close_angles_naive(
        &self,
        vectors: Vec<[f32; 3]>,
        max_angle_rad: f32,
    ) -> PyResult<Vec<([u32; 2], f32)>> {
        let now = Instant::now();
        let res = star::look_up_close_angles_naive(&vectors, max_angle_rad);
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

    let centroids_np = create_2d_numpy_array(py, &centroids);
    let intensities_np = intensities.as_slice().to_pyarray_bound(py);

    Ok((centroids_np, intensities_np))
}

fn create_2d_numpy_array<'py, const N: usize, T>(
    py: Python<'py>,
    points: &[[T; N]],
) -> Bound<'py, numpy::PyArray2<T>>
where
    T: numpy::Element,
{
    let len = points.len(); // n
    let ptr = points.as_ptr() as *const T;
    let total_len = len * N;
    let contiguous_slice = unsafe { std::slice::from_raw_parts(ptr, total_len) };
    let flat_np = contiguous_slice.to_pyarray_bound(py);
    flat_np.reshape((len, N)).unwrap() // Save to unwrap as we know the shape is correct
}

fn numpy_to_vec_3_32f<'py, const L: usize>(
    vectors: &'py numpy::PyReadonlyArray2<'py, f32>,
) -> PyResult<&'py [[f32; L]]> {
    if !vectors.is_c_contiguous() || vectors.ndim() != 2 || vectors.shape()[1] != L {
        return Err(PyRuntimeError::new_err(format!(
            "vectors must be a c_contiguous array with shape=[n, {}]",
            L
        )));
    }
    Ok(util::as_vec_of_arrays(vectors.as_slice()?).unwrap())
}
