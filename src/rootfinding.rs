use super::horner::{horner, horner_eval_f};
use super::{Matrix2, Vector2};

type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;

const PI: f64 = std::f64::consts::PI;

/// The below code defines a struct named Options with three fields: max_iters, tol, and tol_ind.
///
/// Properties:
///
/// * `max_iters`: The `max_iters` property represents the maximum number of iterations allowed for a
/// certain algorithm or process. It is of type `usize`, which means it can only hold non-negative
/// integer values.
/// * `tol`: The `tol` property is a floating-point number that represents the tolerance for convergence
/// in an algorithm. It is used to determine when the algorithm has reached a satisfactory solution.
/// * `tol_ind`: The `tol_ind` property in the `Options` struct represents the tolerance for individual
/// values. It is a floating-point number (`f64`) that determines the acceptable difference between the
/// expected value and the actual value for each element in a calculation or comparison.
#[derive(Debug)]
pub struct Options {
    pub max_iters: usize,
    pub tol: f64,
    pub tol_ind: f64,
}

/// The below code is implementing the `Default` trait for the `Options` struct in Rust. The `Default`
/// trait provides a default value for a type, which can be used when creating an instance of the type
/// without specifying any values. In this case, the `default` function is defined to return an instance
/// of the `Options` struct with default values for the `max_iters`, `tol`, and `tol_ind` fields.
impl Default for Options {
    fn default() -> Self {
        Options {
            max_iters: 2000,
            tol: 1e-12,
            tol_ind: 1e-15,
        }
    }
}

/// The function `make_adjoint` calculates the adjoint matrix between two vectors.
///
/// Arguments:
///
/// * `vr`: A vector representing the direction of the reference frame's x-axis.
/// * `vp`: The parameter `vp` represents a vector `vp = (p, s)`, where `p` and `s` are the components
/// of the vector.
///
/// Returns:
///
/// The function `make_adjoint` returns a `Mat2` object.
#[inline]
pub fn make_adjoint(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    Mat2::new(
        Vector2::<f64>::new(s, -p),
        Vector2::<f64>::new(-p * q, p * r + s),
    )
}

/// The function `make_inverse` calculates the inverse of a 2x2 matrix.
///
/// Arguments:
///
/// * `vr`: A vector representing the row of a 2x2 matrix. The components of the vector are vr.x_ and
/// vr.y_.
/// * `vp`: The parameter `vp` represents a 2D vector with components `x` and `y`.
///
/// Returns:
///
/// The function `make_inverse` returns a `Mat2` object.
#[inline]
pub fn make_inverse(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    let m_adjoint = Mat2::new(
        Vector2::<f64>::new(s, -p),
        Vector2::<f64>::new(-p * q, p * r + s),
    );
    m_adjoint / m_adjoint.det()
}

/// The `delta` function calculates the delta value for the Bairstow's method
///
/// Arguments:
///
/// * `vA`: A vector representing the coefficients of a polynomial equation.
/// * `vr`: The parameter `vr` represents the vector `[-2.0, 0.0]`.
/// * `vp`: The parameter `vp` represents the vector vr - vrj
///
/// Returns:
///
/// The function `delta` returns a `Vec2` object.
///
/// r * p - m   -p
/// q * p       -m
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::delta;
/// use bairstow::vector2::Vector2;
///
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
/// let vd = delta(&vA1, &vri, &vrj);
/// assert_eq!(vd, Vector2::new(0.2, 0.4));
/// ```
#[inline]
pub fn delta(vA: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let mp = make_adjoint(vr, vp); // 2 mul's
    mp.mdot(vA) / mp.det() // 6 mul's + 2 div's
}

/// delta 1 for ri - rj
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::delta1;
/// use bairstow::vector2::Vector2;
///
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, -0.0);
/// let vrj = Vector2::new(4.0, -5.0);
/// let vd = delta1(&vA1, &vri, &vrj);
/// assert_eq!(vd, Vector2::new(0.2, 0.4));
/// ```
#[inline]
pub fn delta1(vA: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    let mp = Matrix2::new(Vec2::new(-s, -p), Vec2::new(p * q, p * r - s));
    mp.mdot(vA) / mp.det() // 6 mul's + 2 div's
}

/// The `suppress_old` function performs zero suppression on a set of vectors.
///
/// Arguments:
///
/// * `vA`: A mutable reference to a Vector2 object representing the coefficients of a polynomial. The
/// coefficients are stored in the x_ and y_ fields of the Vector2 object.
/// * `vA1`: vA1 is a mutable reference to a Vector2 object.
/// * `vri`: The parameter `vri` represents a vector with components `r` and `i`. It is used in the
/// `suppress_old` function to perform calculations.
/// * `vrj`: The parameter `vrj` represents a vector with components `x` and `y`.
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::delta;
/// use bairstow::rootfinding::suppress_old;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut vA = Vector2::new(3.0, 3.0);
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
///
/// suppress_old(&mut vA, &mut vA1, &vri, &vrj);
/// let dr = delta(&vA, &vri, &vA1);
/// assert_approx_eq!(dr.x_, -16.780821917808325);
/// assert_approx_eq!(dr.y_, 1.4383561643835612);
#[inline]
pub fn suppress_old(vA: &mut Vec2, vA1: &mut Vec2, vri: &Vec2, vrj: &Vec2) {
    let (A, B) = (vA.x_, vA.y_);
    let (A1, B1) = (vA1.x_, vA1.y_);
    let vp = vri - vrj;
    let (r, q) = (vri.x_, vri.y_);
    let (p, s) = (vp.x_, vp.y_);
    let f = (r * p) + s;
    let qp = q * p;
    let e = (f * s) - (qp * p);
    let a = ((A * s) - (B * p)) / e;
    let b = ((B * f) - (A * qp)) / e;
    let c = A1 - a;
    let d = (B1 - b) - (a * p);
    vA.x_ = a;
    vA.y_ = b;
    vA1.x_ = ((c * s) - (d * p)) / e;
    vA1.y_ = ((d * f) - (c * qp)) / e;
}

/// The `suppress` function in Rust performs zero suppression on a set of vectors.
///
/// Arguments:
///
/// * `vA`: A vector representing the coefficients of a polynomial function.
/// * `vA1`: The parameter `vA1` is a `Vector2` object representing a vector with two components. It is
/// used as an input parameter in the `suppress` function.
/// * `vri`: The parameter `vri` represents the vector `ri`, and `vrj` represents the vector `rj`. These
/// vectors are used in the calculation of the suppression step in the Bairstow's method for root
/// finding.
/// * `vrj`: The parameter `vrj` represents a vector with coordinates (4.0, 5.0).
/// Zero suppression
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::delta;
/// use bairstow::rootfinding::suppress;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut vA = Vector2::new(3.0, 3.0);
/// let mut vA1 = Vector2::new(1.0, 2.0);
/// let vri = Vector2::new(-2.0, 0.0);
/// let vrj = Vector2::new(4.0, 5.0);
///
/// (vA, vA1) = suppress(&mut vA, &mut vA1, &vri, &vrj);
/// let dr = delta(&vA, &vri, &vA1);
/// assert_approx_eq!(dr.x_, -16.780821917808325);
/// assert_approx_eq!(dr.y_, 1.4383561643835612);
#[inline]
pub fn suppress(vA: &Vec2, vA1: &Vec2, vri: &Vec2, vrj: &Vec2) -> (Vec2, Vec2) {
    let vp = vri - vrj;
    let m_inverse = make_inverse(vri, &vp);
    let va = m_inverse.mdot(vA);
    let mut vc = vA1 - va;
    vc.y_ -= va.x_ * vp.x_;
    let va1 = m_inverse.mdot(&vc);
    (va, va1)
}

/// The `initial_guess` function in Rust calculates the initial guesses for the roots of a polynomial
/// using Bairstow's method.
///
/// Arguments:
///
/// * `coeffs`: A vector of coefficients representing a polynomial.
///
/// Returns:
///
/// The function `initial_guess` returns a vector of `Vector2` structs, which represent the initial
/// guesses for the roots of a polynomial equation.
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::initial_guess;
/// use bairstow::vector2::Vector2;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_guess(&coeffs);
/// ```
pub fn initial_guess(coeffs: &[f64]) -> Vec<Vec2> {
    let mut degree = coeffs.len() - 1;
    let center = -coeffs[1] / (coeffs[0] * degree as f64);
    let centroid = horner_eval_f(coeffs, center); // ???
    let re = centroid.abs().powf(1.0 / (degree as f64));
    degree /= 2;
    degree *= 2; // make even
    let k = PI / (degree as f64);
    let m = center * center + re * re;
    let mut vr0s = Vec::<Vec2>::new();
    for i in (1..degree).step_by(2) {
        let temp = re * (k * i as f64).cos();
        let r0 = 2.0 * (center + temp);
        let t0 = m + 2.0 * center * temp;
        vr0s.push(Vector2::<f64>::new(r0, -t0));
    }
    vr0s
}

/// Parallel Bairstow's method (even degree only)
///
/// The `pbairstow_even` function implements the parallel Bairstow's method for finding roots of
/// even-degree polynomials.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a polynomial.
/// It is assumed that the polynomial has an even degree.
/// * `vrs`: A vector of initial guesses for the roots of the polynomial. Each element of the vector is
/// a complex number representing a root guess.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
/// following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_guess, pbairstow_even, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&coeffs);
/// let (niter, _found) = pbairstow_even(&coeffs, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn pbairstow_even(coeffs: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    let m_rs = vrs.len();
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let mut vri = vrs[i];
            if let Some(tol_i) = pbairstow_even_job(coeffs, i, &mut vri, &mut converged[i], &vrs) {
                if tol < tol_i {
                    tol = tol_i;
                }
            }
            vrs[i] = vri;
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Multi-threading Bairstow's method (even degree only)
///
/// The `pbairstow_even_mt` function implements the multi-threading parallel Bairstow's
/// method for finding roots of even-degree polynomials.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a polynomial.
/// It is assumed that the polynomial has an even degree.
/// * `vrs`: A vector of initial guesses for the roots of the polynomial. Each element of the vector is
/// a complex number representing a root guess.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
/// following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_guess, pbairstow_even_mt, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&coeffs);
/// let (niter, _found) = pbairstow_even_mt(&coeffs, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 8);
/// ```
pub fn pbairstow_even_mt(coeffs: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    use rayon::prelude::*;

    let m_rs = vrs.len();
    let mut vrsc = vec![Vec2::default(); m_rs];
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        vrsc.copy_from_slice(vrs);

        let tol_i = vrs
            .par_iter_mut()
            .zip(converged.par_iter_mut())
            .enumerate()
            .filter(|(_, (_, converged))| !**converged)
            .filter_map(|(i, (vri, converged))| {
                pbairstow_even_job(coeffs, i, vri, converged, &vrsc)
            })
            .reduce(|| tol, |x, y| x.max(y));
        if tol < tol_i {
            tol = tol_i;
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

fn pbairstow_even_job(
    coeffs: &[f64],
    i: usize,
    vri: &mut Vec2,
    converged: &mut bool,
    vrsc: &[Vec2],
) -> Option<f64> {
    let mut pb = coeffs.to_owned();
    // let mut pb = coeffs.to_owned();
    let degree = pb.len() - 1; // degree, assume even
    let mut vA = horner(&mut pb, degree, vri);
    let tol_i = vA.norm_inf();
    if tol_i < 1e-15 {
        *converged = true;
        return None;
    }
    let mut vA1 = horner(&mut pb, degree - 2, vri);
    for (_, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
        // vA1 -= delta(&vA, vrj, &(*vri - vrj));
        suppress_old(&mut vA, &mut vA1, vri, vrj);
    }
    let dt = delta(&vA, vri, &vA1); // Gauss-Seidel fashion
    *vri -= dt;
    Some(tol_i)
}

/// The `initial_autocorr` function calculates the initial guesses for Bairstow's method for finding
/// roots of a polynomial, specifically for the auto-correlation function.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
/// polynomial. The coefficients are ordered from highest degree to lowest degree.
///
/// Returns:
///
/// The function `initial_autocorr` returns a vector of `Vec2` structs.
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::initial_autocorr;
/// use bairstow::vector2::Vector2;
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_autocorr(&coeffs);
/// ```
pub fn initial_autocorr(coeffs: &[f64]) -> Vec<Vec2> {
    let mut degree = coeffs.len() - 1;
    let re = coeffs[degree].abs().powf(1.0 / (degree as f64));
    degree /= 2;
    let k = PI / (degree as f64);
    let m = re * re;
    let mut vr0s = Vec::<Vec2>::new();
    for i in (1..degree).step_by(2) {
        vr0s.push(Vector2::<f64>::new(2.0 * re * (k * i as f64).cos(), -m));
    }
    vr0s
}

/// The `pbairstow_autocorr` function implements the simultaneous Bairstow's method for finding roots of
/// a polynomial, specifically for the auto-correlation function.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
/// polynomial. These coefficients are used to calculate the auto-correlation function.
/// * `vrs`: `vrs` is a vector of complex numbers representing the initial guesses for the roots of the
/// polynomial. Each element of `vrs` is a `Vec2` struct, which contains two fields: `x_` and `y_`.
/// These fields represent the real and imaginary parts of the
/// * `options`: The `Options` struct is used to specify the parameters for the Bairstow's method
/// algorithm. It has the following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_autocorr, pbairstow_autocorr, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&coeffs);
/// let (niter, _found) = pbairstow_autocorr(&coeffs, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 1);
/// ```
pub fn pbairstow_autocorr(coeffs: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    let m_rs = vrs.len();
    let mut converged = vec![false; m_rs];

    for niter in 0..options.max_iters {
        let mut tol = 0.0;

        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let mut vri = vrs[i];
            let tol_i = pbairstow_autocorr_mt_job(coeffs, i, &mut vri, &mut converged[i], &vrs);
            if let Some(tol_i) = tol_i {
                if tol < tol_i {
                    tol = tol_i;
                }
            }
            vrs[i] = vri;
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// The `pbairstow_autocorr_mt` function is a multi-threaded implementation of Bairstow's method for
/// finding roots of a polynomial, specifically for auto-correlation functions.
///
/// Arguments:
///
/// * `coeffs`: The `coeffs` parameter is a slice of `f64` values representing the coefficients of a
/// polynomial. These coefficients are used as input for the Bairstow's method algorithm.
/// * `vrs`: `vrs` is a vector of complex numbers representing the initial guesses for the roots of the
/// polynomial. Each element of `vrs` is a `Vec2` struct, which contains the real and imaginary parts of
/// the complex number.
/// * `options`: The `options` parameter is an instance of the `Options` struct, which contains the
/// following fields:
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_autocorr, pbairstow_autocorr_mt, Options};
///
/// let coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&coeffs);
/// let (niter, _found) = pbairstow_autocorr_mt(&coeffs, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 2);
/// ```
pub fn pbairstow_autocorr_mt(
    coeffs: &[f64],
    vrs: &mut Vec<Vec2>,
    options: &Options,
) -> (usize, bool) {
    use rayon::prelude::*;

    let m_rs = vrs.len();
    let mut vrsc = vec![Vec2::default(); m_rs];
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        vrsc.copy_from_slice(vrs);

        let tol_i = vrs
            .par_iter_mut()
            .zip(converged.par_iter_mut())
            .enumerate()
            .filter(|(_, (_, converged))| !**converged)
            .filter_map(|(i, (vri, converged))| {
                pbairstow_autocorr_mt_job(coeffs, i, vri, converged, &vrsc)
            })
            .reduce(|| tol, |x, y| x.max(y));
        if tol < tol_i {
            tol = tol_i;
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

fn pbairstow_autocorr_mt_job(
    coeffs: &[f64],
    i: usize,
    vri: &mut Vec2,
    converged: &mut bool,
    vrsc: &[Vec2],
) -> Option<f64> {
    let mut pb = coeffs.to_owned();
    // let mut pb = coeffs.to_owned();
    let degree = pb.len() - 1; // assumed divided by 4
    let mut vA = horner(&mut pb, degree, vri);
    let tol_i = vA.norm_inf();
    if tol_i < 1e-15 {
        *converged = true;
        return None;
    }
    let mut vA1 = horner(&mut pb, degree - 2, vri);
    for (_j, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
        // vA1 -= delta(&vA, vrj, &(*vri - vrj));
        suppress_old(&mut vA, &mut vA1, vri, vrj);
        let vrjn = Vector2::<f64>::new(-vrj.x_, 1.0) / vrj.y_;
        // vA1 -= delta(&vA, &vrjn, &(*vri - vrjn));
        suppress_old(&mut vA, &mut vA1, vri, &vrjn);
    }
    let vrin = Vector2::<f64>::new(-vri.x_, 1.0) / vri.y_;
    // vA1 -= delta(&vA, &vrin, &(*vri - vrin));
    suppress_old(&mut vA, &mut vA1, vri, &vrin);
    let dt = delta(&vA, vri, &vA1); // Gauss-Seidel fashion
    *vri -= dt;
    Some(tol_i)
}

/// The `extract_autocorr` function extracts the quadratic function where its roots are within a unit
/// circle.
///
/// x^2 - r*x - t or x^2 + (r/t) * x + (-1/t)
/// (x - a1)(x - a2) = x^2 - (a1 + a2) x + a1 * a2
///
/// Arguments:
///
/// * `vr`: A vector containing two values, representing the coefficients of a quadratic function. The
/// first value represents the coefficient of x^2, and the second value represents the coefficient of x.
///
/// Returns:
///
/// The function `extract_autocorr` returns a `Vec2` struct, which contains two elements `x_` and `y_`.
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::extract_autocorr;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let vr = extract_autocorr(Vector2::new(1.0, -4.0));
///
/// assert_approx_eq!(vr.x_, 0.25);
/// assert_approx_eq!(vr.y_, -0.25);
/// ```
#[allow(dead_code)]
pub fn extract_autocorr(vr: Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: q } = vr;
    let hr = r / 2.0;
    let d = hr * hr + q;
    if d < 0.0 {
        // complex conjugate root
        if q < -1.0 {
            return Vector2::<f64>::new(-r, 1.0) / q;
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
        return Vector2::<f64>::new(a1 + a2, -a1 * a2);
    }
    if a2.abs() > 1.0 {
        a2 = 1.0 / a2;
        return Vector2::<f64>::new(a1 + a2, -a1 * a2);
    }
    // else no need to change
    vr
}
