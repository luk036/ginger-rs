#![allow(non_snake_case)]

use super::horner::{horner_eval_c, horner_eval_f};
use super::Options;
use crate::leja_order::leja_order;
use num_complex::Complex;

const TWO_PI: f64 = std::f64::consts::TAU;

pub fn initial_aberth(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let poly_c = horner_eval_f(coeffs, center);
    let radius = Complex::<f64>::new(-poly_c, 0.0).powf(1.0 / degree as f64);
    (0..degree)
        .map(|i| {
            // note! swap x and y to match C++ Complex{circle2_table_y, circle2_table_x}
            let xcoord = crate::tables::circle2_table_y(i);
            let ycoord = crate::tables::circle2_table_x(i);
            center + radius * Complex::<f64>::new(xcoord, ycoord)
        })
        .collect()
}

/// Initial guess for Aberth's method
///
/// The center $$ c $$ and radius $$ R $$ are derived from the coefficients:
///
/// $$ c = -\frac{a_1}{n a_0}, \qquad R = \sqrt\[n\]{-P(c)}, \qquad z_i = c + R \cdot e^{i\theta_i} $$
///
/// where $$ \theta_i = 2\pi v_i $$ with $$ v_i $$ from a van der Corput sequence
/// for low-discrepancy angular spacing.
///
/// The `initial_aberth` function calculates the initial guesses for Aberth's method given a
/// polynomial's coefficients.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
///   polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
///
/// Returns:
///
/// The function `initial_aberth` returns a vector of `Complex<f64>` values, which represent the initial
/// guesses for the roots of a polynomial.
///
/// # Examples:
///
/// ```
/// use ginger::aberth::initial_aberth_orig;
/// use num_complex::Complex;
/// use approx_eq::assert_approx_eq;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let z0s = initial_aberth_orig(&coeffs);
///
/// assert_approx_eq!(z0s[0].re, 0.6116610247366323);
/// assert_approx_eq!(z0s[0].im, 0.6926747514925476);
/// ```
/// Original initial guess for Aberth's method (equi-angular spacing)
///
/// $$ z_i = c + R \cdot \left(\cos\frac{2\pi k}{n} + i\sin\frac{2\pi k}{n}\right) $$
pub fn initial_aberth_orig(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let poly_c = horner_eval_f(coeffs, center);
    let radius = Complex::<f64>::new(-poly_c, 0.0).powf(1.0 / degree as f64);
    let k = TWO_PI / (degree as f64);
    (0..degree)
        .map(|idx| {
            let theta = k * (0.25 + idx as f64);
            center + radius * Complex::<f64>::new(theta.cos(), theta.sin())
        })
        .collect()
}

fn aberth_job(
    coeffs: &[f64],
    i: usize,
    zi: &mut Complex<f64>,
    zsc: &[Complex<f64>],
    coeffs1: &[f64],
) -> f64 {
    let p_eval = horner_eval_c(coeffs, zi);
    let tol_i = p_eval.l1_norm(); // ???
    let mut p1_eval = horner_eval_c(coeffs1, zi);
    for (_, zj) in zsc.iter().enumerate().filter(|t| t.0 != i) {
        p1_eval -= p_eval / (*zi - zj);
    }
    *zi -= p_eval / p1_eval; // Gauss-Seidel fashion
    tol_i
}

/// Aberth's method
///
/// The `aberth` function implements Aberth's method for finding roots of a polynomial.
///
/// Each estimate $$ z_i $$ is updated by the correction:
///
/// $$ z_i' = z_i - \frac{P(z_i)}{P'(z_i)} $$
///
/// where the derivative is modified by all other root estimates:
///
/// $$ P'(z_i) = P_1(z_i) - \sum_{\substack{j=1\\j\neq i}}^n \frac{P(z_i)}{z_i - z_j} $$
///
/// Here $$ P_1(z) $$ is the ordinary derivative of $$ P(z) $$. Convergence is cubic.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
///   polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
/// * `zs`: A vector of complex numbers representing the initial guesses for the roots of the polynomial.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::Options;
/// use ginger::aberth::{initial_aberth, aberth};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&coeffs);
/// let (niter, _found) = aberth(&coeffs, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn aberth(coeffs: &[f64], zs: &mut [Complex<f64>], options: &Options) -> (usize, bool) {
    let m_zs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let coeffs1: Vec<_> = coeffs[0..degree]
        .iter()
        .enumerate()
        .map(|(i, ci)| ci * (degree - i) as f64)
        .collect();

    for niter in 0..options.max_iters {
        let mut tolerance = 0.0;

        for i in 0..m_zs {
            let mut zi = zs[i];
            let tol_i = aberth_job(coeffs, i, &mut zi, zs, &coeffs1);
            if tolerance < tol_i {
                tolerance = tol_i;
            }
            zs[i] = zi;
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Multi-threading Aberth's method
///
/// The `aberth_mt` function in Rust implements the multi-threaded Aberth's method for root finding.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///   polynomial. The polynomial is defined by the equation:
/// * `zs`: A mutable reference to a vector of Complex numbers. These numbers represent the initial
///   guesses for the roots of the polynomial equation.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///   following fields:
///
/// # Examples:
///
/// ```
/// use ginger::rootfinding::Options;
/// use ginger::aberth::{initial_aberth, aberth_mt};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&coeffs);
/// let (niter, _found) = aberth_mt(&coeffs, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 6);
/// ```
pub fn aberth_mt(coeffs: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    fn aberth_job2(
        coeffs: &[f64],
        i: usize,
        zi: &mut Complex<f64>,
        zsc: &[Complex<f64>],
        coeffs1: &[f64],
    ) -> f64 {
        let p_eval = horner_eval_c(coeffs, zi);
        let tol_i = p_eval.l1_norm(); // ???
        let mut p1_eval = horner_eval_c(coeffs1, zi);
        for (_, zj) in zsc.iter().enumerate().filter(|t| t.0 != i) {
            p1_eval -= p_eval / (*zi - zj);
        }
        *zi -= p_eval / p1_eval; // Gauss-Seidel fashion
        tol_i
    }

    use rayon::prelude::*;
    let m_zs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let coeffs1: Vec<_> = (0..degree)
        .map(|i| coeffs[i] * (degree - i) as f64)
        .collect();
    let mut zsc = vec![Complex::default(); m_zs];

    for niter in 0..options.max_iters {
        let mut tolerance = 0.0;
        zsc.copy_from_slice(zs);

        let tol_i = zs
            .par_iter_mut()
            .enumerate()
            .map(|(i, zi)| aberth_job2(coeffs, i, zi, &zsc, &coeffs1))
            .reduce(|| tolerance, |x, y| x.max(y));
        if tolerance < tol_i {
            tolerance = tol_i;
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Initial guess for Aberth's method using auto-correlation
///
/// $$ R = \sqrt\[n\]{|a_n|},\qquad R \leftarrow \max(R, 1/R),\qquad z_i = c + R \cdot e^{i\theta_i} $$
///
/// The `initial_aberth_autocorr` function calculates initial guesses for Aberth's method
/// specifically tailored for auto-correlation polynomials.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients
///   of a polynomial. The coefficients are ordered from highest degree to lowest degree.
///
/// Returns:
///
/// The function returns a vector of `Complex<f64>` values, representing the initial guesses
/// for the roots of a polynomial.
pub fn initial_aberth_autocorr(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1; // assume even
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let poly_c = horner_eval_f(coeffs, center);
    let mut radius = poly_c.abs().powf(1.0 / degree as f64);
    if radius > 1.0 {
        radius = 1.0 / radius;
    }
    (0..degree / 2)
        .map(|i| {
            // note! swap x and y to match C++ Complex{circle2_table_y, circle2_table_x}
            let xcoord = crate::tables::circle2_table_y(i);
            let ycoord = crate::tables::circle2_table_x(i);
            center + radius * Complex::<f64>::new(xcoord, ycoord)
        })
        .collect()
}

/// Aberth's method job for auto-correlation polynomials
///
/// This internal function performs a single iteration of Aberth's method for auto-correlation
/// polynomials, considering both the root and its reciprocal.
///
/// Arguments:
///
/// * `coeffs`: Polynomial coefficients
/// * `i`: Current root index
/// * `zi`: Current root value (mutable)
/// * `zsc`: Current approximations of all roots
/// * `coeffs1`: Derivative coefficients
///
/// Returns:
///
/// The tolerance value for convergence checking.
fn aberth_autocorr_job(
    coeffs: &[f64],
    i: usize,
    zi: &mut Complex<f64>,
    zsc: &[Complex<f64>],
    coeffs1: &[f64],
) -> f64 {
    let p_eval = horner_eval_c(coeffs, zi);
    let tol_i = p_eval.l1_norm(); // ???
    let mut p1_eval = horner_eval_c(coeffs1, zi);
    for (_, zj) in zsc.iter().enumerate().filter(|t| t.0 != i) {
        p1_eval -= p_eval / (*zi - zj);
        p1_eval -= p_eval / (*zi - 1.0 / zj);
    }
    *zi -= p_eval / p1_eval; // Gauss-Seidel fashion
    tol_i
}

/// Aberth's method for auto-correlation polynomials
///
/// Each estimate $$ z_i $$ is updated by the correction, considering both $$ z $$
/// and $$ 1/\bar{z} $$ in the derivative:
///
/// $$ P'(z_i) = P_1(z_i) - \sum_{\substack{j=1\\j\neq i}}^n
///    \left(\frac{P(z_i)}{z_i - z_j} + \frac{P(z_i)}{z_i - 1/z_j}\right) $$
///
/// The `aberth_autocorr` function implements Aberth's method specifically for
/// auto-correlation polynomials, where roots come in reciprocal pairs.
///
/// Arguments:
///
/// * `coeffs`: The polynomial coefficients (highest to lowest degree)
/// * `zs`: Mutable slice of complex root approximations
/// * `options`: Iteration options (max iterations, tolerance)
///
/// Returns:
///
/// A tuple of (number of iterations, whether convergence was achieved)
pub fn aberth_autocorr(
    coeffs: &[f64],
    zs: &mut [Complex<f64>],
    options: &Options,
) -> (usize, bool) {
    let m_zs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let coeffs1: Vec<_> = coeffs[0..degree]
        .iter()
        .enumerate()
        .map(|(i, ci)| ci * (degree - i) as f64)
        .collect();

    for niter in 0..options.max_iters {
        let mut tolerance = 0.0;

        for i in 0..m_zs {
            let mut zi = zs[i];
            let tol_i = aberth_autocorr_job(coeffs, i, &mut zi, zs, &coeffs1);
            if tolerance < tol_i {
                tolerance = tol_i;
            }
            zs[i] = zi;
        }
        if tolerance < options.tolerance {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Reconstruct a monic polynomial from its roots using Leja ordering
///
/// $$ P(x) = \prod_{i=1}^n (x - r_i) = x^n + a_{n-1} x^{n-1} + \cdots + a_0 $$
///
/// The coefficients are computed by repeated convolution. Starting from $$ \[1\] $$,
/// for each root $$ r $$:
///
/// $$ c_{k+1} \leftarrow c_{k+1} - r \cdot c_k $$
///
/// Given a set of complex roots, reconstruct the monic polynomial coefficients
/// (highest degree first) by multiplying (x - root) factors. Leja ordering is
/// applied for numerical accuracy.
///
/// Arguments:
///
/// * `zs` - Input vector of complex roots
///
/// Returns:
///
/// Monic polynomial coefficients (highest degree first)
pub fn poly_from_roots(zs: &[Complex<f64>]) -> Vec<f64> {
    if zs.is_empty() {
        return vec![1.0];
    }
    let ordered = leja_order(zs.to_vec());
    let mut coeffs = vec![Complex::new(1.0, 0.0)];
    for z in &ordered {
        let mut prev = coeffs[0];
        for item in coeffs.iter_mut().skip(1) {
            let old = *item;
            *item -= z * prev;
            prev = old;
        }
        coeffs.push(-z * prev);
    }
    coeffs.iter().map(|c| c.re).collect()
}

/// Reconstruct a monic polynomial from its autocorrelation roots
///
/// Auto-correlation (palindromic) polynomials have roots in reciprocal pairs.
/// The aberth_autocorr functions find the degree/2 "independent" roots.
/// This function adds the reciprocal of each root (1/z) to get the full set
/// of degree roots, then reconstructs with Leja ordering.
///
/// Arguments:
///
/// * `zs` - Roots found by aberth_autocorr
///
/// Returns:
///
/// Monic polynomial coefficients (highest degree first)
pub fn poly_from_autocorr_roots(zs: &[Complex<f64>]) -> Vec<f64> {
    if zs.is_empty() {
        return vec![1.0];
    }
    // Add reciprocals to account for the palindromic root-pair structure
    let mut all_roots: Vec<Complex<f64>> = Vec::with_capacity(2 * zs.len());
    for z in zs {
        all_roots.push(*z);
        all_roots.push(1.0 / z);
    }
    poly_from_roots(&all_roots)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_horner_eval() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let z = Complex::new(0.0, 0.0);
        let p_eval = horner_eval_c(&coeffs, &z);
        assert_eq!(p_eval.re, 10.0);
        assert_eq!(p_eval.im, 0.0);
        let z = Complex::new(1.0, 0.0);
        let p_eval = horner_eval_c(&coeffs, &z);
        assert_eq!(p_eval.re, 576.0);
        assert_eq!(p_eval.im, 0.0);
    }

    #[test]
    fn test_aberth() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut zrs = initial_aberth(&coeffs);
        let (niter, found) = aberth(&coeffs, &mut zrs, &Options::default());
        assert_eq!(niter, 5);
        assert!(found);
    }

    #[test]
    fn test_aberth_mt() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut zrs = initial_aberth(&coeffs);
        let (niter, found) = aberth_mt(&coeffs, &mut zrs, &Options::default());
        assert_eq!(niter, 6);
        assert!(found);
    }

    #[test]
    fn test_aberth_autocorr() {
        let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
        let mut zrs = initial_aberth_autocorr(&coeffs);
        let (niter, found) = aberth_autocorr(&coeffs, &mut zrs, &Options::default());
        assert!(niter <= 7);
        assert!(found);
    }

    #[test]
    fn test_poly_from_roots() {
        // Polynomial (x-1)(x-2) = x^2 - 3x + 2
        let roots = vec![Complex::new(1.0, 0.0), Complex::new(2.0, 0.0)];
        let coeffs = poly_from_roots(&roots);
        assert_eq!(coeffs.len(), 3);
        assert!((coeffs[0] - 1.0).abs() < 1e-12);
        assert!((coeffs[1] + 3.0).abs() < 1e-12);
        assert!((coeffs[2] - 2.0).abs() < 1e-12);
    }

    #[test]
    fn test_poly_from_autocorr_roots() {
        // Simple test: just check it doesn't panic and returns right length
        let roots = vec![Complex::new(0.5, 0.5), Complex::new(0.5, -0.5)];
        let coeffs = poly_from_autocorr_roots(&roots);
        // Should have 2*roots.len() + 1 coefficients
        assert_eq!(coeffs.len(), 5);
        // Monic
        assert!((coeffs[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_poly_from_roots_empty() {
        let coeffs = poly_from_roots(&[]);
        assert_eq!(coeffs, vec![1.0]);
    }
}
