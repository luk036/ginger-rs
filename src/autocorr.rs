//! Auto-correlation polynomial root-finding (palindromic polynomials).
//!
//! Mirrors `autocorr.hpp` / `autocorr.cpp` in ginger-cpp: the palindromic
//! variant of Bairstow's method where roots appear in reciprocal pairs.

use crate::execution_policy::{atomic_decoupled_run, jacobi_mt_run, sequential_run, Step};
use crate::rootfinding::{delta, horner, roots_from_quadratic, suppress_old, Options, Vec2};
use crate::seqlock::AtomicVec2;
use num_complex::Complex;

/// One Bairstow Newton correction respecting palindromic (auto-correlation)
/// symmetry: each neighbor factor contributes both `vrj` and its reciprocal
/// image `(-vrj.x, 1) / vrj.y`; the factor itself also suppresses its own
/// reciprocal image. Shared by the sequential, Jacobi-MT and atomic execution
/// policies.
///
/// This mirrors `ginger::detail::autocorr_bairstow_step` in ginger-cpp.
pub struct AutocorrBairstowStep<'a> {
    /// Polynomial coefficients (highest degree first).
    pub coeffs: &'a [f64],
    /// Degree of the polynomial.
    pub degree: usize,
    /// Convergence options.
    pub options: &'a Options,
}

impl Step<Vec2> for AutocorrBairstowStep<'_> {
    type Cell = AtomicVec2;

    fn run<G, N>(&self, idx: usize, get: G, neighbors: N) -> (f64, Option<Vec2>)
    where
        G: Fn(usize) -> Vec2,
        N: Iterator<Item = usize>,
    {
        let vri = get(idx);
        let mut coeffs1 = self.coeffs.to_owned(); // horner corrupts the array
        let mut v_big_a = horner(&mut coeffs1, self.degree, &vri);
        let tol_i = v_big_a.norm_inf();
        if tol_i < self.options.tol_ind {
            return (0.0, None);
        }
        let mut v_big_a1 = horner(&mut coeffs1, self.degree - 2, &vri);
        for j in neighbors {
            let vrj = get(j);
            suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrj);
            let vrjn = Vec2::new(-vrj.x_, 1.0) / vrj.y_;
            suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrjn);
        }
        let vrin = Vec2::new(-vri.x_, 1.0) / vri.y_;
        suppress_old(&mut v_big_a, &mut v_big_a1, &vri, &vrin);
        let dt = delta(&v_big_a, &vri, &v_big_a1); // Gauss-Seidel fashion
        (tol_i, Some(vri - dt))
    }
}

/// The `initial_autocorr` function calculates the initial guesses for Bairstow's method for finding
/// roots of a polynomial, specifically for the auto-correlation function.
///
/// $$ R = \sqrt\[n\]{|a_n|},\qquad R \leftarrow \max(R, 1/R),\qquad m = n/2 $$
/// $$ \theta_k = \frac{k\pi}{m},\qquad (r_k, q_k) = (2R\cos\theta_k,\; -R^2) $$
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. The coefficients are ordered from highest degree to lowest degree.
///
/// Returns:
///
/// The function `initial_autocorr` returns a vector of `Vec2` structs.
///
/// # Examples:
///
/// ```
/// use ginger::autocorr::initial_autocorr;
/// use ginger::vector2::Vector2;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_autocorr(&coeffs);
/// ```
pub fn initial_autocorr(coeffs: &[f64]) -> Vec<Vec2> {
    let degree = coeffs.len() - 1;
    let radius = coeffs[degree].abs().powf(1.0 / (degree as f64));
    let degree = degree / 2;
    let m = radius * radius;
    let num_points = degree / 2;
    (0..num_points)
        .map(|i| Vec2::new(2.0 * radius * crate::tables::cos_pi_vdc2(i), -m))
        .collect()
}

/// The `pbairstow_autocorr` function implements the simultaneous Bairstow's method for finding roots of
/// a polynomial, specifically for the auto-correlation function.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. These coefficients are used to calculate the auto-correlation function.
/// * `vrs`: `vrs` is a vector of complex numbers representing the initial guesses for the roots of the
///   polynomial. Each element of `vrs` is a `Vec2` struct, which contains two fields: `x_` and `y_`.
///   These fields represent the real and imaginary parts of the
/// * `options`: The `Options` struct is used to specify the parameters for the Bairstow's method
///   algorithm. It has the following fields:
///
/// # Examples:
///
/// ```
/// use ginger::autocorr::{initial_autocorr, pbairstow_autocorr};
/// use ginger::rootfinding::Options;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&coeffs);
/// let (niter, found) = pbairstow_autocorr(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_autocorr(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    let step = AutocorrBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    sequential_run(vrs, options, &step)
}

/// The `pbairstow_autocorr_mt` function is a multi-threaded implementation of Bairstow's method for
/// finding roots of a polynomial, specifically for auto-correlation functions.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. These coefficients are used as input for the Bairstow's method algorithm.
/// * `vrs`: `vrs` is a vector of complex numbers representing the initial guesses for the roots of the
///   polynomial. Each element of `vrs` is a `Vec2` struct, which contains the real and imaginary parts of
///   the complex number.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::autocorr::{initial_autocorr, pbairstow_autocorr_mt};
/// use ginger::rootfinding::Options;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&coeffs);
/// let (niter, found) = pbairstow_autocorr_mt(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_autocorr_mt(coeffs: &[f64], vrs: &mut [Vec2], options: &Options) -> (usize, bool) {
    let step = AutocorrBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    jacobi_mt_run(vrs, options, &step)
}

/// Atomic multi-threading Bairstow's method (auto-correlation, decoupled)
///
/// The `pbairstow_autocorr_atomic` function is a multi-threaded implementation of Bairstow's
/// method for finding roots of a polynomial with auto-correlation (palindromic) symmetry using a
/// single atomic working buffer.
///
/// Unlike `pbairstow_autocorr_mt` (Jacobi snapshot + per-iteration barrier), the atomic buffer is
/// built once: each thread owns a chunk of factor slots (single-writer, multi-reader via a
/// seqlock) and runs its own iteration loop independently. There is no per-iteration
/// synchronization — a thread exits as soon as its own chunk converges or the maximum number of
/// iterations is exceeded. The iteration count is therefore NON-DETERMINISTIC (the returned count
/// is the maximum across threads), and `found` is true only if every chunk converged.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. These coefficients are used as input for the Bairstow's method algorithm.
/// * `vrs`: `vrs` is a vector of complex numbers representing the initial guesses for the roots of
///   the polynomial. Each element of `vrs` is a `Vec2` struct, which contains the real and
///   imaginary parts of the complex number.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::autocorr::{initial_autocorr, pbairstow_autocorr_atomic};
/// use ginger::rootfinding::Options;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&coeffs);
/// let (niter, found) = pbairstow_autocorr_atomic(&coeffs, &mut vrs, &Options::default());
///
/// assert!(niter > 0);
/// assert!(found);
/// ```
pub fn pbairstow_autocorr_atomic(
    coeffs: &[f64],
    vrs: &mut [Vec2],
    options: &Options,
) -> (usize, bool) {
    let step = AutocorrBairstowStep {
        coeffs,
        degree: coeffs.len() - 1,
        options,
    };
    atomic_decoupled_run(vrs, options, &step)
}

/// The `extract_autocorr` function extracts quadratic factors from a polynomial with auto-correlation
/// property.
///
/// Given a quadratic $$ x^2 - r x - q $$, computes its roots and replaces
/// any root $$ |z| > 1 $$ with its reciprocal $$ 1/z $$. The normalized
/// factor is recovered from the adjusted roots via Vieta:
///
/// $$ r' = z_1' + z_2', \qquad q' = -z_1' z_2' $$
///
/// where $$ z_k' = z_k $$ if $$ |z_k| \le 1 $$, else $$ z_k' = 1/z_k $$.
///
/// Arguments:
///
/// * `vr`: A vector containing two values, representing the coefficients of a quadratic function. The
///   first value represents the coefficient of x^2, and the second value represents the coefficient of x.
///
/// Returns:
///
/// The function `extract_autocorr` returns a `Vec2` struct, which contains two elements `x_` and `y_`.
///
/// # Examples:
///
/// ```
/// use ginger::autocorr::extract_autocorr;
/// use ginger::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let vr = extract_autocorr(Vector2::new(1.0, -4.0));
///
/// assert_approx_eq!(vr.x_, 0.25);
/// assert_approx_eq!(vr.y_, -0.25);
/// ```
pub fn extract_autocorr(vr: Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: q } = vr;
    let hr = r / 2.0;
    let d = hr * hr + q;
    if d < 0.0 {
        // complex conjugate root
        if q < -1.0 {
            return Vec2::new(-r, 1.0) / q;
        }
    }
    // two real roots
    let mut a1 = hr + (if hr >= 0.0 { d.sqrt() } else { -d.sqrt() });
    let mut a2 = -q / a1;

    if a1.abs() > 1.0 {
        if a2.abs() > 1.0 {
            a2 = 1.0 / a2;
        }
        a1 = 1.0 / a1;
        return Vec2::new(a1 + a2, -a1 * a2);
    }
    if a2.abs() > 1.0 {
        a2 = 1.0 / a2;
        return Vec2::new(a1 + a2, -a1 * a2);
    }
    // else no need to change
    vr
}

/// Reconstruct a monic polynomial from its autocorrelation quadratic factors
///
/// Auto-correlation (palindromic) polynomials have roots in reciprocal pairs.
/// Each quadratic factor $x^2 - r x - q$ found by `pbairstow_autocorr` carries 2 roots.
/// This function adds the reciprocal of each root, then reconstructs the full
/// monic polynomial with Leja ordering for numerical accuracy.
///
/// $$ P(x) = \prod_{i=1}^{m} (x^2 - r_i x - q_i)(x^{-2} - r_i x^{-1} - q_i) $$
///
/// Arguments:
///
/// * `vrs` - Quadratic factors from pbairstow_autocorr
///
/// Returns:
///
/// Monic polynomial coefficients (highest degree first)
pub fn poly_from_autocorr_factors(vrs: &[Vec2]) -> Vec<f64> {
    if vrs.is_empty() {
        return vec![1.0];
    }
    // Each factor x^2 - r*x - q contributes 2 roots. For palindromic/autocorrelation
    // polynomials, the reciprocal of each root is also a root. Collect all roots
    // and their reciprocals, then reconstruct with Leja ordering.
    let mut all_roots: Vec<Complex<f64>> = Vec::with_capacity(4 * vrs.len());
    for vr in vrs {
        let (r1, r2) = roots_from_quadratic(vr);
        all_roots.push(r1);
        all_roots.push(r2);
        all_roots.push(1.0 / r1);
        all_roots.push(1.0 / r2);
    }
    crate::aberth::poly_from_roots(&all_roots)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx_eq::assert_approx_eq;

    #[test]
    fn test_initial_autocorr() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let guesses = initial_autocorr(&coeffs);

        assert_eq!(guesses.len(), 2);
        // Verify the first guess is reasonable
        assert!(guesses[0].x_.abs() > 0.0);
        assert!(guesses[0].y_.abs() > 0.0);
    }

    #[test]
    fn test_pbairstow_autocorr() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_autocorr(&coeffs);
        let options = Options::default();

        let (niter, found) = pbairstow_autocorr(&coeffs, &mut vrs, &options);

        assert!(niter > 0);
        assert!(found);
        // Verify at least one root is close to actual root
        let mut has_root = false;
        for vr in vrs {
            let val = horner(&mut coeffs.clone(), coeffs.len() - 1, &vr);
            if val.norm_inf() < options.tolerance {
                has_root = true;
                break;
            }
        }
        assert!(has_root);
    }

    #[test]
    fn test_pbairstow_autocorr_atomic() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_autocorr(&coeffs);
        let (niter, found) = pbairstow_autocorr_atomic(&coeffs, &mut vrs, &Options::default());
        assert!(niter > 0);
        assert!(found);
    }

    const FIR_COEFFS: [f64; 49] = [
        -0.00196191,
        -0.00094597,
        -0.00023823,
        0.00134667,
        0.00380494,
        0.00681596,
        0.0097864,
        0.01186197,
        0.0121238,
        0.00985211,
        0.00474894,
        -0.00281751,
        -0.01173923,
        -0.0201885,
        -0.02590168,
        -0.02658216,
        -0.02035729,
        -0.00628271,
        0.01534627,
        0.04279982,
        0.0732094,
        0.10275561,
        0.12753013,
        0.14399228,
        0.15265722,
        0.14399228,
        0.12753013,
        0.10275561,
        0.0732094,
        0.04279982,
        0.01534627,
        -0.00628271,
        -0.02035729,
        -0.02658216,
        -0.02590168,
        -0.0201885,
        -0.01173923,
        -0.00281751,
        0.00474894,
        0.00985211,
        0.0121238,
        0.01186197,
        0.0097864,
        0.00681596,
        0.00380494,
        0.00134667,
        -0.00023823,
        -0.00094597,
        -0.00196191,
    ];

    #[test]
    fn test_pbairstow_autocorr_atomic_fir() {
        // Decoupled threads under load (libtest runs tests in parallel) can need
        // many iterations to converge; give a generous iteration budget.
        let options = Options {
            max_iters: 20000,
            tolerance: 1e-2,
            ..Options::default()
        };
        let mut vrs = initial_autocorr(&FIR_COEFFS);
        let (_, found) = pbairstow_autocorr_atomic(&FIR_COEFFS, &mut vrs, &options);
        assert!(found);
    }

    #[test]
    fn test_pbairstow_autocorr_atomic_reconstruction() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut vrs = initial_autocorr(&coeffs);
        let (_, found) = pbairstow_autocorr_atomic(&coeffs, &mut vrs, &Options::default());
        assert!(found);
        let monic = poly_from_autocorr_factors(&vrs);
        let scale = coeffs[0];
        for (i, c) in coeffs.iter().enumerate() {
            assert!(
                (monic[i] * scale - c).abs() < 1e-8,
                "coefficient {i} mismatch: {} vs {}",
                monic[i] * scale,
                c
            );
        }
    }

    #[test]
    fn test_extract_autocorr() {
        let vr = Vec2::new(1.0, -4.0);
        let result = extract_autocorr(vr);

        assert_approx_eq!(result.x_, 0.25);
        assert_approx_eq!(result.y_, -0.25);
    }

    #[test]
    fn test_poly_from_autocorr_factors() {
        let vrs = vec![Vec2::new(3.0, -2.0)];
        let coeffs = poly_from_autocorr_factors(&vrs);
        // With reciprocals: roots are 2, 1, 0.5, 1.0
        // polynomial = (x-2)(x-1)(x-0.5)(x-1) = ...
        assert_eq!(coeffs.len(), 5);
        assert!((coeffs[0] - 1.0).abs() < 1e-12);
    }
}
