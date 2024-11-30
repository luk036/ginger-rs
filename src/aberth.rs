#![allow(non_snake_case)]

use super::Options;
use lds_rs::lds::Circle;
use num_complex::Complex;

const TWO_PI: f64 = std::f64::consts::TAU;

/// Horner evalution (float)
///
/// The `horner_eval_f` function in Rust implements the Horner's method for evaluating a polynomial with
/// given coefficients at a specific value.
///
/// Arguments:
///
/// * `coeffs`: A vector of floating-point coefficients representing a polynomial. The coefficients are
///             ordered from highest degree to lowest degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 +
///             94x^5 + 150x^4 + 94x^
/// * `zval`: The `zval` parameter in the `horner_eval_f` function represents the value at which the
///             polynomial is evaluated. It is of type `f64`, which means it is a floating-point number.
///
/// Returns:
///
/// The function `horner_eval_f` returns a `f64` value, which is the result of evaluating the polynomial
/// with the given coefficients at the specified value `zval`.
///
/// # Examples:
///
/// ```
/// use ginger::aberth::horner_eval_f;
/// use approx_eq::assert_approx_eq;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_f(&coeffs, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// ```
pub fn horner_eval_f(coeffs: &[f64], zval: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, coeff| acc * zval + coeff)
}

/// Horner evalution (complex)
///
/// The `horner_eval_c` function in Rust implements the Horner evaluation method for complex
/// polynomials.
///
/// Arguments:
///
/// * `coeffs`: A vector of coefficients representing a polynomial. The coefficients are in descending
///             order of degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 + 94x^5 + 150x^4 + 94x^3 + 75
/// * `zval`: The `zval` parameter is a complex number that represents the value at which the polynomial is evaluated.
///
/// Returns:
///
/// The function `horner_eval_c` returns a complex number of type `Complex<f64>`.
///
/// # Examples:
///
/// ```
/// use ginger::aberth::horner_eval_c;
/// use approx_eq::assert_approx_eq;
/// use num_complex::Complex;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_c(&coeffs, &Complex::new(1.0, 2.0));
///
/// assert_approx_eq!(px.re, 6080.0);
/// assert_approx_eq!(px.im, 9120.0);
/// ```
pub fn horner_eval_c(coeffs: &[f64], zval: &Complex<f64>) -> Complex<f64> {
    coeffs
        .iter()
        .fold(Complex::<f64>::new(0.0, 0.0), |acc, coeff| {
            acc * zval + coeff
        })
}

pub fn initial_aberth(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let poly_c = horner_eval_f(coeffs, center);
    let radius = Complex::<f64>::new(-poly_c, 0.0).powf(1.0 / degree as f64);
    let mut c_gen = Circle::new(2);
    (0..degree)
        .map(|_idx| {
            let [ycoord, xcoord] = c_gen.pop(); // Note: y, x
            center + radius * Complex::<f64>::new(xcoord, ycoord)
        })
        .collect()
}

/// Initial guess for Aberth's method
///
/// The `initial_aberth` function calculates the initial guesses for Aberth's method given a
/// polynomial's coefficients.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///             polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
///             polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
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
/// <pre>
///                 P ⎛z ⎞
///      new          ⎝ i⎠
///     z    = z  - ───────
///      i      i   P' ⎛z ⎞
///                    ⎝ i⎠
/// where
///                           degree
///                         _____
///                         ╲
///                          ╲    P ⎛z ⎞
///                           ╲     ⎝ i⎠
///     P' ⎛z ⎞ = P  ⎛z ⎞ -   ╱   ───────
///        ⎝ i⎠    1 ⎝ i⎠    ╱    z  - z
///                         ╱      i    j
///                         ‾‾‾‾‾
///                         j ≠ i
/// </pre>
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
///             polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
///             polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
/// * `zs`: A vector of complex numbers representing the initial guesses for the roots of the polynomial.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///             following fields:
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
///             polynomial. The polynomial is defined by the equation:
/// * `zs`: A mutable reference to a vector of Complex numbers. These numbers represent the initial
///             guesses for the roots of the polynomial equation.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
///             following fields:
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

pub fn initial_aberth_autocorr(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1; // assume even
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let poly_c = horner_eval_f(coeffs, center);
    let mut radius = poly_c.abs().powf(1.0 / degree as f64);
    if radius > 1.0 {
        radius = 1.0 / radius;
    }
    let mut c_gen = Circle::new(2);
    (0..degree / 2)
        .map(|_idx| {
            let [y, x] = c_gen.pop();
            center + radius * Complex::<f64>::new(x, y)
        })
        .collect()
}

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
}
