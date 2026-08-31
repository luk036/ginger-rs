//! Execution policies for parallelizable root-finding solvers.
//!
//! Strategy + Template-Method decomposition:
//!
//! - Each algorithm supplies a per-root *step* (the `Step` trait) that computes one
//!   Newton correction for a root/factor. The step reads the current values of
//!   the other roots through a `get` closure provided by the policy, and
//!   returns the corrected value; the policy owns all writes to its data
//!   layout (live state, frozen snapshot, or seqlock buffer).
//! - Each execution policy owns the iteration loop, the scheduling of the
//!   steps, and the convergence aggregation.
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

/// A per-root Newton correction step shared by every execution policy.
///
/// `run(idx, get, neighbors)` computes one correction for root/factor `idx`:
///
/// - `get(j)` reads the current value of root/factor `j` from the policy's
///   data layout (live state, frozen snapshot, or seqlock buffer).
/// - `neighbors` enumerates the indices of the other roots/factors to
///   suppress (ascending order, or round-robin for the atomic variant).
///
/// Returns `(tol, Some(new_value))` when the root was updated, or
/// `(0.0, None)` when it is already converged (below `options.tol_ind`).
///
/// This mirrors the per-algorithm step functors (`even_bairstow_step`,
/// `autocorr_bairstow_step`, `aberth_step`, ...) in ginger-cpp's
/// `execution_policy.hpp`.
pub trait Step<T>: Sync {
    /// Atomic cell type used by the decoupled atomic policy for this step.
    type Cell: AtomicCell<T>;

    /// Compute one Newton correction for root/factor `idx`.
    fn run<G, N>(&self, idx: usize, get: G, neighbors: N) -> (f64, Option<T>)
    where
        G: Fn(usize) -> T,
        N: Iterator<Item = usize>;
}

/// Sequential Gauss-Seidel execution: roots are updated in-place in ascending
/// index order within each iteration (later roots see earlier updates).
///
/// Converged roots contribute a zero tolerance and are not written; their
/// values stay frozen as neighbors for the remaining sweeps.
pub fn sequential_run<T, S>(state: &mut [T], options: &Options, step: &S) -> (usize, bool)
where
    T: Copy,
    S: Step<T>,
{
    let m = state.len();
    for niter in 0..options.max_iters {
        let mut tolerance: f64 = 0.0;
        for idx in 0..m {
            let (tol_i, new_value) = step.run(idx, |j| state[j], (0..m).filter(move |&j| j != idx));
            tolerance = tolerance.max(tol_i);
            if let Some(value) = new_value {
                state[idx] = value;
            }
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Jacobi multi-threaded execution: each iteration reads a frozen snapshot of
/// the roots and updates them in parallel, so the iteration is
/// order-independent. Small problems (<= 4 roots) fall back to the sequential
/// policy.
pub fn jacobi_mt_run<T, S>(state: &mut [T], options: &Options, step: &S) -> (usize, bool)
where
    T: Copy + Default + Send + Sync,
    S: Step<T> + Sync,
{
    if !should_parallelize(state.len()) {
        return sequential_run(state, options, step);
    }
    let m = state.len();
    let mut snapshot = vec![T::default(); m];
    for niter in 0..options.max_iters {
        snapshot.copy_from_slice(state);

        let updates: Vec<(f64, Option<T>)> = (0..m)
            .into_par_iter()
            .map(|idx| step.run(idx, |j| snapshot[j], (0..m).filter(move |&j| j != idx)))
            .collect();

        let mut tolerance: f64 = 0.0;
        for (idx, (tol_i, new_value)) in updates.into_iter().enumerate() {
            tolerance = tolerance.max(tol_i);
            if let Some(value) = new_value {
                state[idx] = value;
            }
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Atomic decoupled execution: a single seqlock buffer is built once; each
/// thread owns a chunk of slots (single-writer, multi-reader) and runs its own
/// iteration loop independently with no per-iteration barrier, reading the
/// other slots in round-robin order to reduce contention. The returned
/// iteration count is the maximum across threads and is non-deterministic.
pub fn atomic_decoupled_run<T, S>(state: &mut [T], options: &Options, step: &S) -> (usize, bool)
where
    T: Copy + Send + Sync,
    S: Step<T>,
{
    let num_roots = state.len();
    let buffer: Vec<S::Cell> = state.iter().copied().map(S::Cell::new).collect();
    let buffer_ref = &buffer;

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
                        // Round-robin suppression order: each thread reads the
                        // other slots in a different rotation, reducing
                        // concurrent access to the same slot.
                        let neighbors = (0..num_roots - 1).map(move |k| (idx + k + 1) % num_roots);
                        let (tol_i, new_value) = step.run(idx, |j| buffer_ref[j].load(), neighbors);
                        max_tol = max_tol.max(tol_i);
                        if let Some(value) = new_value {
                            buffer_ref[idx].store(value);
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
