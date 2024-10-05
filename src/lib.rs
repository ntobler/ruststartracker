use numpy::{self, PyUntypedArrayMethods};
use pyo3::{
    exceptions::PyRuntimeError, pyclass, pyfunction, pymethods, pymodule, types::PyModule,
    wrap_pyfunction, Bound, PyRef, PyRefMut, PyResult,
};
use std::{time::Instant, usize};

mod star;
mod tree;
mod trianglefinder;
mod util;

#[pymodule]
fn libruststartracker(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(main_scratch, m)?)?;
    m.add_class::<TriangleFinder>()?;
    m.add_class::<IterTriangleFinder>()?;
    m.add_class::<StarMatcher>()?;
    m.add_class::<UnitVectorLookup>()?;
    Ok(())
}

#[pyfunction]
fn main_scratch() -> PyResult<()> {
    // TODO remove me
    let connections_ab = vec![[234, 5643], [1, 2], [2, 4], [3, 9], [2, 6]];
    let connections_ac = vec![[345, 2343], [8, 2], [3, 4], [1, 7], [0, 5], [3, 1]];
    let connections_bc = vec![[435, 4355], [1, 0], [4, 8], [8, 1], [1, 9]];

    let finder =
        trianglefinder::TriangleFinder::new(connections_ab, connections_ac, connections_bc);

    let now = Instant::now();

    let found: Option<[u32; 3]> = finder.get();
    println!("Time passed: {:?}", now.elapsed());
    match found {
        None => println!("not found"),
        Some(value) => println!("found {:?}", value),
    }

    println!("And now we iterate");

    let iter_finder = trianglefinder::IterTriangleFinder::new(finder);

    for value in iter_finder {
        println!("found {:?}", value)
    }

    let a: Vec<[f32; 3]> = vec![
        [-0.36230828, -0.88712497, -0.28590558],
        [-0.00418702, -0.98547485, -0.16976978],
        [-0.3743498, -0.87157188, -0.31658281],
    ];
    let b: Vec<[f32; 3]> = vec![
        [0.20531651, -0.0681248, 0.97632174],
        [-0.16902335, 0.03930809, 0.98482789],
        [0.2206221, -0.10030539, 0.97018798],
    ];

    println!("---");
    let now2 = Instant::now();
    let svd = star::attitude_svd(&a, &b);
    println!("Time passed: {:?}\n", now2.elapsed());

    match svd {
        None => println!("not found"),
        Some(value) => println!("found {:?}", value),
    }

    let vs = vec![
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.0, 0.0, -1.0],
    ];
    let l = tree::UnitVectorLookup::new(&vs);
    let res = l.lookup_nearest(&[0.0, 0.0, -1.0]);
    print!("look up result is {:?}", res);

    println!("what");

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
    ) -> Self {
        StarMatcher {
            inner: star::StarMatcher::new(
                stars_xyz,
                max_inter_star_angle,
                inter_star_angle_tolerance,
                n_minimum_matches,
            ),
        }
    }

    pub fn find(
        &self,
        obs_xyz: Vec<[f32; 3]>,
    ) -> PyResult<([f32; 4], Vec<u32>, u32, Vec<[f32; 3]>, f32)> {
        let now = Instant::now();
        let res = self.inner.find(obs_xyz);
        let duration_s = now.elapsed().as_secs_f32();
        match res {
            None => Err(PyRuntimeError::new_err("Could not find attitude.")),
            Some(x) => Ok((x.quat, x.match_ids, x.n_matches, x.obs_matched, duration_s)),
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
    ) -> PyResult<(Vec<[u32; 2]>, Vec<f32>, Vec<f32>)> {
        let now = Instant::now();
        let res = star::get_inter_star_index(&self.inner, &vectors, max_angle_rad);
        println!("Time passed: {:?}", now.elapsed());
        Ok(res)
    }

    pub fn get_inter_star_index_numpy<'py>(
        &self,
        vectors: numpy::PyReadonlyArray2<'py, f32>,
        max_angle_rad: f32,
    ) -> PyResult<(Vec<[u32; 2]>, Vec<f32>, Vec<f32>)> {
        let now = Instant::now();
        let vectors_inner = numpy_to_vec_3_32f(&vectors).unwrap();
        let res = star::get_inter_star_index(&self.inner, vectors_inner, max_angle_rad);
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

fn numpy_to_vec_3_32f<'py, const L: usize>(
    vectors: &'py numpy::PyReadonlyArray2<'py, f32>,
) -> PyResult<&[[f32; L]]> {
    if !vectors.is_c_contiguous() || vectors.ndim() != 2 || vectors.shape()[1] != L {
        return Err(PyRuntimeError::new_err(format!(
            "vectors must be a c_contiguous array with shape=[n, {}]",
            L
        )));
    }
    Ok(util::as_vec_of_arrays(vectors.as_slice()?).unwrap())
}
