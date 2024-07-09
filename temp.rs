#![allow(non_snake_case)]

use super::Options;
use num_complex::Complex;

const TWO_PI: f64 = std::f64::consts::TAU;

pub fn horner_eval_f(coeffs: &[f64], zval: f64) -> f64 {
    coeffs.iter().fold(0.0, |acc, coeff| acc * zval + coeff)
}

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
    let k = TWO_PI / (degree as f64);
    (0..degree)
        .map(|idx| {
            let theta = k * (0.25 + idx as f64);
            center + radius * Complex::<f64>::new(theta.cos(), theta.sin())
        })
        .collect()
}

/// Aberth's method
pub fn aberth(coeffs: &[f64], zs: &mut [Complex<f64>], options: &Options) -> (usize, bool) {
    let m_zs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
                                   // let coeffs1: Vec<_> = (0..degree)
                                   //     .map(|i| coeffs[i] * (degree - i) as f64)
                                   //     .collect();
    let coeffs1: Vec<_> = coeffs[0..degree]
        .iter()
        .enumerate()
        .map(|(i, ci)| ci * (degree - i) as f64)
        .collect();
    let mut converged = vec![false; m_zs];

    for niter in 0..options.max_iters {
        let mut tolerance = 0.0;

        for i in 0..m_zs {
            if converged[i] {
                continue;
            }
            let mut zi = zs[i];
            if let Some(tol_i) = aberth_job(coeffs, i, &mut zi, &mut converged[i], zs, &coeffs1) {
                if tolerance < tol_i {
                    tolerance = tol_i;
                }
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
pub fn aberth_mt(coeffs: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    use rayon::prelude::*;

    let m_zs = zs.len();
    let degree = coeffs.len() - 1; // degree, assume even
    let coeffs1: Vec<_> = (0..degree)
        .map(|i| coeffs[i] * (degree - i) as f64)
        .collect();
    let mut zsc = vec![Complex::default(); m_zs];
    let mut converged = vec![false; m_zs];

    for niter in 0..options.max_iters {
        let mut tolerance = 0.0;
        zsc.copy_from_slice(zs);

        let tol_i = zs
            .par_iter_mut()
            .zip(converged.par_iter_mut())
            .enumerate()
            .filter(|(_, (_, converged))| !**converged)
            .filter_map(|(i, (zi, converged))| aberth_job(coeffs, i, zi, converged, &zsc, &coeffs1))
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

fn aberth_job(
    coeffs: &[f64],
    i: usize,
    zi: &mut Complex<f64>,
    converged: &mut bool,
    zsc: &[Complex<f64>],
    coeffs1: &[f64],
) -> Option<f64> {
    let p_eval = horner_eval_c(coeffs, zi);
    let tol_i = p_eval.l1_norm(); // ???
    if tol_i < 1e-15 {
        *converged = true;
        return None;
    }
    let mut p1_eval = horner_eval_c(coeffs1, zi);
    for (_, zj) in zsc.iter().enumerate().filter(|t| t.0 != i) {
        p1_eval -= p_eval / (*zi - zj);
    }
    *zi -= p_eval / p1_eval; // Gauss-Seidel fashion
    Some(tol_i)
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
}
