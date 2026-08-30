//! Facade entry points that auto-select the execution policy.
//!
//! Mirrors `solve.hpp` in ginger-cpp and `solve.py` in ginger (Python): each
//! solver exposes a facade that dispatches to the single-threaded,
//! multi-threaded, or atomic variant based on the problem size
//! ([`should_parallelize`]) or an explicit [`SolveMode`].

use crate::aberth::{aberth, aberth_atomic, aberth_autocorr, aberth_mt};
use crate::execution_policy::should_parallelize;
use crate::rootfinding::{
    pbairstow_autocorr, pbairstow_autocorr_atomic, pbairstow_autocorr_mt, pbairstow_even,
    pbairstow_even_atomic, pbairstow_even_mt, Options,
};
use num_complex::Complex;

/// Execution policy selector for the solver facades.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolveMode {
    /// Pick multi-threaded when [`should_parallelize`] holds, else sequential.
    Automatic,
    Sequential,
    MultiThreaded,
    Atomic,
}

/// Aberth-Ehrlich solver with automatic policy selection.
///
/// `AUTOMATIC` dispatches to the multi-threaded variant when
/// [`should_parallelize`] holds, otherwise to the single-threaded variant.
pub fn solve_aberth(
    coeffs: &[f64],
    zs: &mut Vec<Complex<f64>>,
    options: &Options,
    mode: SolveMode,
) -> (usize, bool) {
    match mode {
        SolveMode::Sequential => aberth(coeffs, zs, options),
        SolveMode::MultiThreaded => aberth_mt(coeffs, zs, options),
        SolveMode::Atomic => aberth_atomic(coeffs, zs, options),
        SolveMode::Automatic => {
            if should_parallelize(zs.len()) {
                aberth_mt(coeffs, zs, options)
            } else {
                aberth(coeffs, zs, options)
            }
        }
    }
}

/// Bairstow solver (even degree) with automatic policy selection.
pub fn solve_pbairstow_even(
    coeffs: &[f64],
    vrs: &mut Vec<crate::rootfinding::Vec2>,
    options: &Options,
    mode: SolveMode,
) -> (usize, bool) {
    match mode {
        SolveMode::Sequential => pbairstow_even(coeffs, vrs, options),
        SolveMode::MultiThreaded => pbairstow_even_mt(coeffs, vrs, options),
        SolveMode::Atomic => pbairstow_even_atomic(coeffs, vrs, options),
        SolveMode::Automatic => {
            if should_parallelize(vrs.len()) {
                pbairstow_even_mt(coeffs, vrs, options)
            } else {
                pbairstow_even(coeffs, vrs, options)
            }
        }
    }
}

/// Bairstow solver for autocorrelation polynomials with policy selection.
pub fn solve_pbairstow_autocorr(
    coeffs: &[f64],
    vrs: &mut Vec<crate::rootfinding::Vec2>,
    options: &Options,
    mode: SolveMode,
) -> (usize, bool) {
    match mode {
        SolveMode::Sequential => pbairstow_autocorr(coeffs, vrs, options),
        SolveMode::MultiThreaded => pbairstow_autocorr_mt(coeffs, vrs, options),
        SolveMode::Atomic => pbairstow_autocorr_atomic(coeffs, vrs, options),
        SolveMode::Automatic => {
            if should_parallelize(vrs.len()) {
                pbairstow_autocorr_mt(coeffs, vrs, options)
            } else {
                pbairstow_autocorr(coeffs, vrs, options)
            }
        }
    }
}

/// Aberth solver for autocorrelation polynomials.
///
/// Only the single-threaded variant exists in this crate, so every mode
/// delegates to [`aberth_autocorr`].
pub fn solve_aberth_autocorr(
    coeffs: &[f64],
    zs: &mut [Complex<f64>],
    options: &Options,
    _mode: SolveMode,
) -> (usize, bool) {
    aberth_autocorr(coeffs, zs, options)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_should_parallelize_threshold() {
        assert!(!should_parallelize(0));
        assert!(!should_parallelize(4));
        assert!(should_parallelize(5));
    }

    #[test]
    fn test_solve_aberth_automatic_small() {
        // 4 roots -> automatic dispatches to the single-threaded variant.
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut zrs = crate::aberth::initial_aberth(&coeffs);
        let (_, found) = solve_aberth(&coeffs, &mut zrs, &Options::default(), SolveMode::Automatic);
        assert!(found);
    }

    #[test]
    fn test_solve_aberth_explicit_modes() {
        let coeffs = vec![
            1.0, 0.0, 2.0, 0.0, 3.0, 0.0, 4.0, 0.0, 3.0, 0.0, 2.0, 0.0, 1.0,
        ];
        let mut zrs = crate::aberth::initial_aberth(&coeffs);
        let (_, found) = solve_aberth(
            &coeffs,
            &mut zrs,
            &Options::default(),
            SolveMode::MultiThreaded,
        );
        assert!(found);
        let mut zrs = crate::aberth::initial_aberth(&coeffs);
        let (_, found) = solve_aberth(&coeffs, &mut zrs, &Options::default(), SolveMode::Atomic);
        assert!(found);
    }

    #[test]
    fn test_solve_pbairstow_even() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = crate::rootfinding::initial_guess(&coeffs);
        let (_, found) =
            solve_pbairstow_even(&coeffs, &mut vrs, &Options::default(), SolveMode::Automatic);
        assert!(found);
    }
}
