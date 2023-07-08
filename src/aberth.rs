#![allow(non_snake_case)]

use super::Options;
use num::Complex;
// use lds_rs::lds::Circle;

const TWO_PI: f64 = std::f64::consts::TAU;

/// Horner evalution (float)
/// 
/// The `horner_eval_f` function in Rust implements the Horner's method for evaluating a polynomial with
/// given coefficients at a specific value.
/// 
/// Arguments:
/// 
/// * `coeffs`: A vector of floating-point coefficients representing a polynomial. The coefficients are
/// ordered from highest degree to lowest degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 +
/// 94x^5 + 150x^4 + 94x^
/// * `zval`: The `zval` parameter in the `horner_eval_f` function represents the value at which the
/// polynomial is evaluated. It is of type `f64`, which means it is a floating-point number.
/// 
/// Returns:
/// 
/// The function `horner_eval_f` returns a `f64` value, which is the result of evaluating the polynomial
/// with the given coefficients at the specified value `zval`.
///
/// # Examples:
///
/// ```
/// use bairstow::aberth::horner_eval_f;
/// use approx_eq::assert_approx_eq;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_f(&coeffs, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// ```
pub fn horner_eval_f(coeffs: &[f64], zval: f64) -> f64 {
    let mut res = coeffs[0];
    for coeff in coeffs.iter().skip(1) {
        res = res * zval + coeff;
    }
    res
}

/// Horner evalution (complex)
/// 
/// The `horner_eval_c` function in Rust implements the Horner evaluation method for complex
/// polynomials.
/// 
/// Arguments:
/// 
/// * `coeffs`: A vector of coefficients representing a polynomial. The coefficients are in descending
/// order of degree. For example, the polynomial 10x^8 + 34x^7 + 75x^6 + 94x^5 + 150x^4 + 94x^3 + 75
/// * `zval`: The `zval` parameter is a complex number that represents the value at which the polynomial
/// is evaluated.
/// 
/// Returns:
/// 
/// The function `horner_eval_c` returns a complex number of type `Complex<f64>`.
///
/// # Examples:
///
/// ```
/// use bairstow::aberth::horner_eval_c;
/// use approx_eq::assert_approx_eq;
/// use num::Complex;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval_c(&coeffs, &Complex::new(1.0, 2.0));
///
/// assert_approx_eq!(px.re, 6080.0);
/// assert_approx_eq!(px.im, 9120.0);
/// ```
pub fn horner_eval_c(coeffs: &[f64], zval: &Complex<f64>) -> Complex<f64> {
    let mut res = Complex::<f64>::new(coeffs[0], 0.0);
    for coeff in coeffs.iter().skip(1) {
        res *= zval;
        res += coeff;
    }
    res
}

/// Initial guess for Aberth's method
/// 
/// The `initial_aberth` function calculates the initial guesses for Aberth's method given a
/// polynomial's coefficients.
/// 
/// Arguments:
/// 
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
/// polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
/// polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
/// 
/// Returns:
/// 
/// The function `initial_aberth` returns a vector of `Complex<f64>` values, which represent the initial
/// guesses for the roots of a polynomial.
///
/// # Examples:
///
/// ```
/// use bairstow::aberth::initial_aberth;
/// use num::Complex;
/// use approx_eq::assert_approx_eq;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let z0s = initial_aberth(&coeffs);
///
/// assert_approx_eq!(z0s[0].re, 0.6116610247366323);
/// assert_approx_eq!(z0s[0].im, 0.6926747514925476);
/// ```
pub fn initial_aberth(coeffs: &[f64]) -> Vec<Complex<f64>> {
    let degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let Pc = horner_eval_f(coeffs, center);
    let re = Complex::<f64>::new(-Pc, 0.0).powf(1.0 / degree as f64);
    let k = TWO_PI / (degree as f64);
    let mut z0s = vec![];
    for idx in 0..degree {
        let theta = k * (0.25 + idx as f64);
        let z0 = center + re * Complex::<f64>::new(theta.cos(), theta.sin());
        z0s.push(z0);
    }
    z0s
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
/// polynomial. The coefficients are ordered from highest degree to lowest degree. For example, if the
/// polynomial is `3x^2 + 2x + 1`, the `coeffs` slice would
/// * `zs`: A vector of complex numbers representing the initial guesses for the roots of the
/// polynomial.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
/// following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::Options;
/// use bairstow::aberth::{initial_aberth, aberth};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&coeffs);
/// let (niter, _found) = aberth(&coeffs, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn aberth(coeffs: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    let m_rs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let mut converged = vec![false; m_rs];
    let mut pb = vec![0.0; degree];
    for i in 0..degree {
        pb[i] = coeffs[i] * (degree - i) as f64;
    }
    for niter in 0..options.max_iters {
        let mut tol = 0.0;
        let mut rx = vec![];

        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let mut job = || {
                let zi = &zs[i];
                let pp = horner_eval_c(coeffs, zi);
                let tol_i = pp.l1_norm(); // ???
                if tol_i < 1e-15 {
                    converged[i] = true;
                    rx.push(tol_i);
                }
                let mut pp1 = horner_eval_c(&pb, zi);
                for (_, zj) in zs.iter().enumerate().filter(|t| t.0 != i) {
                    pp1 -= pp / (zi - zj);
                }
                zs[i] -= pp / pp1; // Gauss-Seidel fashion
                rx.push(tol_i);
            };
            job();
        }
        for result in rx.iter() {
            if tol < *result {
                tol = *result;
            }
        }
        if tol < options.tol {
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
/// polynomial. The polynomial is defined by the equation:
/// * `zs`: A mutable reference to a vector of Complex numbers. These numbers represent the initial
/// guesses for the roots of the polynomial equation.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
/// following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::Options;
/// use bairstow::aberth::{initial_aberth, aberth_mt};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&coeffs);
/// let (niter, _found) = aberth_mt(&coeffs, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 7);
/// ```
pub fn aberth_mt(coeffs: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use std::sync::Arc;
    use threadpool::ThreadPool;

    let n_workers = 4; // assume 4 cores

    let m_rs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let mut pb = vec![0.0; degree];
    for i in 0..degree {
        pb[i] = coeffs[i] * (degree - i) as f64;
    }
    // let mut zsc = zs.clone();
    let coeffs = coeffs; // make imutatable
    let pa_share = Arc::new(coeffs.to_owned());
    let pb_share = Arc::new(pb.to_owned());
    // let zs_share = Arc::new(Mutex::new(&zs));

    let mut converged = vec![false; m_rs];

    for niter in 0..options.max_iters {
        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);
        let mut n_jobs = 0;

        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let tx = tx.clone();
            let zsc = zs.clone();
            // let pac = coeffs.clone();
            // let zi = Complex::<f64>::default();
            let pa_clone = Arc::clone(&pa_share);
            let pb_clone = Arc::clone(&pb_share);
            // let zs_clone = Arc::clone(&zs_share);

            n_jobs += 1;
            pool.execute(move || {
                let zi = zsc[i];
                let pp = horner_eval_c(&pa_clone, &zi);
                let tol_i = pp.l1_norm(); // ???
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                } else {
                    let mut pp1 = horner_eval_c(&pb_clone, &zi);
                    for (j, zj) in zsc.iter().enumerate() {
                        // exclude i
                        if j == i {
                            continue;
                        }
                        pp1 -= pp / (zi - zj);
                    }
                    let dt = pp / pp1; // Gauss-Seidel fashion
                    tx.send((Some((tol_i, dt)), i))
                        .expect("channel will be there waiting for a pool");
                }
            });
        }
        // let mut zsw = zs_share.lock().unwrap();
        for (res, i) in rx.iter().take(n_jobs) {
            if let Some(result) = res {
                let (toli, dt) = result;
                if tol < toli {
                    tol = toli;
                }
                zs[i] -= dt;
            } else {
                converged[i] = true;
            }
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}
