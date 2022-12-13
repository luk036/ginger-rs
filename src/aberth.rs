#![allow(non_snake_case)]

use super::Options;
use num::Complex;
// use lds_rs::lds::Circle;

const TWO_PI: f64 = std::f64::consts::TAU;

/// Horner evalution (float)
///
/// Examples:
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
/// Examples:
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
/// Examples:
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
    let ppc = horner_eval_f(coeffs, center);
    let re = Complex::<f64>::new(-ppc, 0.0).powf(1.0 / degree as f64);
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
/// <pre>
///                 P ⎛z ⎞
///      new          ⎝ i⎠
///     z    = z  - ───────
///      i      i   P' ⎛z ⎞
///                    ⎝ i⎠
/// where
///                           n
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
/// Examples:
///
/// ```
/// use bairstow::rootfinding::Options;
/// use bairstow::aberth::{initial_aberth, aberth};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&pa);
/// let (niter, _found) = aberth(&pa, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn aberth(pa: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    let m = zs.len();
    let n = pa.len() - 1; // degree, assume even
    let mut converged = vec![false; m];
    let mut pb = vec![0.0; n];
    for i in 0..n {
        pb[i] = pa[i] * (n - i) as f64;
    }
    for niter in 0..options.max_iter {
        let mut tol = 0.0;
        let mut rx = vec![];

        for i in 0..m {
            if converged[i] {
                continue;
            }
            let mut job = || {
                let zi = &zs[i];
                let pp = horner_eval_c(pa, zi);
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
    (options.max_iter, false)
}

/// Multi-threading Aberth's method
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::Options;
/// use bairstow::aberth::{initial_aberth, aberth_th};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut zrs = initial_aberth(&pa);
/// let (niter, _found) = aberth_th(&pa, &mut zrs, &Options::default());
///
/// assert_eq!(niter, 7);
/// ```
pub fn aberth_th(pa: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use std::sync::Arc;
    use threadpool::ThreadPool;

    let n_workers = 4; // assume 4 cores

    let m = zs.len();
    let n = pa.len() - 1; // degree, assume even
    let mut coeffs = vec![0.0; n];
    let n = pa.len() - 1; // degree, assume even
    for k in 0..n {
        coeffs[k] = pa[k] * (n - k) as f64;
    }
    // let mut zsc = zs.clone();
    let coeffs = coeffs; // make imutatable
    let pa_share = Arc::new(pa.to_owned());
    let pb_share = Arc::new(coeffs);
    // let zs_share = Arc::new(Mutex::new(&zs));

    let mut converged = vec![false; m];

    for niter in 0..options.max_iter {
        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);
        let mut n_jobs = 0;

        for i in 0..m {
            if converged[i] {
                continue;
            }
            let tx = tx.clone();
            let zsc = zs.clone();
            // let pac = pa.clone();
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
    (options.max_iter, false)
}
