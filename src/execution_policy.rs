//! Execution policies for parallelizable root-finding solvers.
//!
//! Strategy + Template-Method decomposition:
//!
//! - Each algorithm supplies a per-root *job* closure that computes one Newton
//!   correction and returns the residual tolerance (or `None` when the root is
//!   already converged).
//! - Each execution policy owns the iteration loop, the scheduling of the jobs,
//!   and the convergence aggregation.
//!
//! This mirrors `execution_policy.hpp` in ginger-cpp and `_policy.py` in
//! ginger (Python): the algorithm (Aberth/Bairstow variant) is decoupled from
//! the execution mode (sequential Gauss-Seidel, Jacobi-MT, or atomic
//! decoupled).

use crate::seqlock::{AtomicComplex, AtomicVec2};
use crate::Options;
use num_complex::Complex;
use rayon::prelude::*;
use std::sync::Mutex;

/// Number of roots above which the multi-threaded policies are used by default.
pub const PARALLEL_THRESHOLD: usize = 4;

/// Whether `num_roots` should use the multi-threaded execution policy.
#[inline]
pub fn should_parallelize(num_roots: usize) -> bool {
    num_roots > PARALLEL_THRESHOLD
}

/// Atomic cell abstraction over a pair of `f64` values (seqlock-backed).
///
/// Implemented by [`AtomicComplex`] and [`AtomicVec2`] so the atomic policy can
/// be generic over the state element type.
pub trait AtomicCell<T>: Sync {
    /// Create a cell initialized with `value`.
    fn new(value: T) -> Self;
    /// Torn-free read of the current value.
    fn load(&self) -> T;
    /// Store a new value; single-writer per slot.
    fn store(&self, value: T);
}

impl AtomicCell<Complex<f64>> for AtomicComplex {
    #[inline]
    fn new(value: Complex<f64>) -> Self {
        AtomicComplex::new(value)
    }
    #[inline]
    fn load(&self) -> Complex<f64> {
        self.load()
    }
    #[inline]
    fn store(&self, value: Complex<f64>) {
        self.store(value);
    }
}

impl AtomicCell<Vector2F64> for AtomicVec2 {
    #[inline]
    fn new(value: Vector2F64) -> Self {
        AtomicVec2::new(value)
    }
    #[inline]
    fn load(&self) -> Vector2F64 {
        self.load()
    }
    #[inline]
    fn store(&self, value: Vector2F64) {
        self.store(value);
    }
}

type Vector2F64 = crate::vector2::Vector2<f64>;

/// Sequential Gauss-Seidel execution: roots are updated in-place in ascending
/// index order within each iteration (later roots see earlier updates).
///
/// `job(i, zi, converged, state)` computes the correction for root `i`, reading
/// the live `state` for neighbours. Returns `Some(tol)` when an update was made
/// (contributing to the tolerance), or `None` when the root is already
/// converged (setting `converged`).
///
/// `from` is the iteration count of the first sweep (0- or 1-based, matching
/// the caller's convention).
pub fn sequential_run<T, F>(
    state: &mut [T],
    options: &Options,
    from: usize,
    job: F,
) -> (usize, bool)
where
    T: Copy,
    F: Fn(usize, &mut T, &mut bool, &[T]) -> Option<f64> + Sync,
{
    let m = state.len();
    let mut converged = vec![false; m];
    for niter in from..options.max_iters {
        let mut tolerance: f64 = 0.0;
        for i in 0..m {
            if converged[i] {
                continue;
            }
            let mut zi = state[i];
            if let Some(tol_i) = job(i, &mut zi, &mut converged[i], state) {
                tolerance = tolerance.max(tol_i);
            }
            state[i] = zi;
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Jacobi multi-threaded execution: each iteration reads a frozen snapshot of
/// the roots and updates them in parallel; converged roots are skipped.
///
/// `from` is the iteration count of the first sweep (0- or 1-based).
pub fn jacobi_mt_run<T, F>(state: &mut [T], options: &Options, from: usize, job: F) -> (usize, bool)
where
    T: Copy + Default + Send + Sync,
    F: Fn(usize, &mut T, &mut bool, &[T]) -> Option<f64> + Sync,
{
    let m = state.len();
    let mut snapshot = vec![T::default(); m];
    let mut converged = vec![false; m];
    for niter in from..options.max_iters {
        let mut tolerance: f64 = 0.0;
        snapshot.copy_from_slice(state);

        let tol_i = state
            .par_iter_mut()
            .zip(converged.par_iter_mut())
            .enumerate()
            .filter(|(_, (_, converged))| !**converged)
            .filter_map(|(i, (zi, converged))| job(i, zi, converged, &snapshot))
            .reduce(|| tolerance, |x, y| x.max(y));
        tolerance = tolerance.max(tol_i);
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Atomic decoupled execution: a single seqlock buffer is built once; each
/// thread owns a chunk of slots (single-writer, multi-reader) and runs its own
/// iteration loop independently with no per-iteration barrier. The returned
/// iteration count is the maximum across threads and is non-deterministic.
pub fn atomic_decoupled_run<T, C, F>(state: &mut [T], options: &Options, job: F) -> (usize, bool)
where
    T: Copy + Send + Sync,
    C: AtomicCell<T> + Sync,
    F: Fn(usize, &[C]) -> Option<f64> + Sync,
{
    let num_roots = state.len();
    let buffer: Vec<C> = state.iter().copied().map(C::new).collect();
    let buffer_ref = &buffer;
    let job_ref = &job;

    let use_mt = should_parallelize(num_roots);
    let num_threads = if use_mt {
        rayon::current_num_threads().min(num_roots)
    } else {
        1
    };
    let chunk_size = (num_roots + num_threads - 1) / num_threads;

    let results: Mutex<Vec<(usize, bool)>> = Mutex::new(Vec::new());
    let results_ref = &results;
    rayon::scope(|s| {
        for t in 0..num_threads {
            let start = t * chunk_size;
            let end = (start + chunk_size).min(num_roots);
            if start >= end {
                break;
            }
            s.spawn(move |_| {
                let mut niter = 0;
                loop {
                    if niter == options.max_iters {
                        results_ref.lock().unwrap().push((options.max_iters, false));
                        return;
                    }
                    let mut max_tol: f64 = 0.0;
                    for idx in start..end {
                        if let Some(tol_i) = job_ref(idx, buffer_ref) {
                            max_tol = max_tol.max(tol_i);
                        }
                    }
                    if max_tol < options.tolerance {
                        results_ref.lock().unwrap().push((niter, true));
                        return;
                    }
                    niter += 1;
                }
            });
        }
    });

    let mut niter_max = 0;
    let mut all_converged = true;
    for (niter, converged) in results.into_inner().unwrap() {
        niter_max = niter_max.max(niter);
        all_converged &= converged;
    }

    for (i, z) in state.iter_mut().enumerate() {
        *z = buffer[i].load();
    }
    (niter_max, all_converged)
}
