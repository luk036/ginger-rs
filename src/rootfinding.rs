use super::{Matrix2, Vector2};

type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;

const PI: f64 = std::f64::consts::PI;

#[derive(Debug)]
pub struct Options {
    pub max_iters: usize,
    pub tol: f64,
    pub tol_ind: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            max_iters: 2000,
            tol: 1e-12,
            tol_ind: 1e-15,
        }
    }
}

/// Make adjoint between two vectors
#[inline]
pub fn make_adjoint(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    Mat2::new(
        Vector2::<f64>::new(s, -p),
        Vector2::<f64>::new(-p * q, p * r + s),
    )
}

/**
 * @brief
 *
 * @param vr
 * @param vp
 * @return Mat2
 */
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

/// Extract the quadratic function where its roots are within a unit circle
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

/// Zero suppression (original)
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

/// Horner evalution
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::horner_eval;
/// use approx_eq::assert_approx_eq;
///
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval(&mut coeffs, 8, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// assert_approx_eq!(coeffs[3], 460.0);
/// ```
#[inline]
pub fn horner_eval(coeffs: &mut [f64], degree: usize, zval: f64) -> f64 {
    for idx in 0..degree {
        coeffs[idx + 1] += coeffs[idx] * zval;
    }
    coeffs[degree]
}

/// Horner evalution for Bairstow's method
///
/// # Examples:
///
/// ```
/// use bairstow::rootfinding::horner;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut coeffs = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner(&mut coeffs, 8, &Vector2::new(-1.0, -2.0));
///
/// assert_approx_eq!(px.x_, 114.0);
/// assert_approx_eq!(px.y_, 134.0);
/// assert_approx_eq!(coeffs[3], 15.0);           
/// ```
pub fn horner(coeffs: &mut [f64], degree: usize, vr: &Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: q } = vr;
    for idx in 0..(degree - 1) {
        coeffs[idx + 1] += coeffs[idx] * r;
        coeffs[idx + 2] += coeffs[idx] * q;
    }
    Vector2::<f64>::new(coeffs[degree - 1], coeffs[degree])
}

/// Initial guess for Bairstow's method
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
    let mut pb = coeffs.to_owned();
    let centroid = horner_eval(&mut pb, degree, center); // ???
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
    let degree = coeffs.len() - 1; // degree, assume even
    let m_rs = vrs.len();
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        for i in 0..m_rs {
            if converged[i] {
                continue;
            }
            let mut pb = coeffs.to_owned();
            let vri = vrs[i];
            let vA = &mut horner(&mut pb, degree, &vri);
            let tol_i = vA.norm_inf();
            if tol_i < 1e-15 {
                converged[i] = true;
                continue;
            } else {
                let vA1 = &mut horner(&mut pb, degree - 2, &vri);
                if tol < tol_i {
                    tol = tol_i;
                }
                for (j, vrj) in vrs.iter().enumerate().take(m_rs) {
                    if j == i {
                        continue;
                    }
                    // vA1 -= delta(vA, vrj, &(vri - vrj));
                    suppress_old(vA, vA1, &vri, vrj);
                }
                vrs[i] -= delta(vA, &vri, vA1); // Gauss-Seidel fashion
            }
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
    use std::sync::mpsc::channel;
    use threadpool::ThreadPool;

    let n_workers = 4; // assume 4 cores
    let m_rs = vrs.len();
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);
        let mut n_jobs = 0;
        for i in (0..m_rs).filter(|x| !converged[*x]) {
            // if converged[i] {
            //     continue;
            // }
            let tx = tx.clone();
            let vrsc = vrs.clone();
            let mut pb = coeffs.to_owned();

            n_jobs += 1;
            pool.execute(move || {
                // let mut pb = coeffs.to_owned();
                let degree = pb.len() - 1; // degree, assume even
                let vri = vrsc[i];
                let mut vA = horner(&mut pb, degree, &vri);
                let tol_i = vA.norm_inf();
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                } else {
                    let mut vA1 = horner(&mut pb, degree - 2, &vri);
                    for (_, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
                        // vA1 -= delta(&vA, vrj, &(vri - vrj));
                        suppress_old(&mut vA, &mut vA1, &vri, vrj);
                    }
                    let dt = delta(&vA, &vri, &vA1); // Gauss-Seidel fashion
                    tx.send((Some((tol_i, dt)), i))
                        .expect("channel will be there waiting for a pool");
                }
            });
        }
        for (res, i) in rx.iter().take(n_jobs) {
            if let Some(result) = res {
                let (toli, dt) = result;
                if tol < toli {
                    tol = toli;
                }
                vrs[i] -= dt;
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

/// Initial guess for Bairstow's method (specific for auto-correlation function)
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

/// Simultenous Bairstow's method (specific for auto-correlation function)
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
            let mut job = || {
                let mut pb = coeffs.to_owned();
                let degree = coeffs.len() - 1; // assumed divided by 4
                let vri = vrs[i];
                let mut vA = horner(&mut pb, degree, &vri);
                let tol_i = vA.norm_inf();
                if tol_i < 1e-15 {
                    converged[i] = true;
                    return tol_i;
                }
                let mut vA1 = horner(&mut pb, degree - 2, &vri);
                for (_j, vrj) in vrs.iter().enumerate().filter(|t| t.0 != i) {
                    // vA1 -= delta(&vA, vrj, &(vri - vrj));
                    suppress_old(&mut vA, &mut vA1, &vri, vrj);
                    let vrjn = Vector2::<f64>::new(-vrj.x_, 1.0) / vrj.y_;
                    // vA1 -= delta(&vA, &vrjn, &(vri - vrjn));
                    suppress_old(&mut vA, &mut vA1, &vri, &vrjn);
                }
                let vrin = Vector2::<f64>::new(-vri.x_, 1.0) / vri.y_;
                // vA1 -= delta(&vA, &vrin, &(vri - vrin));
                suppress_old(&mut vA, &mut vA1, &vri, &vrin);
                vrs[i] -= delta(&vA, &vri, &vA1); // Gauss-Seidel fashion
                tol_i
            };
            let tol_i = job();
            if tol < tol_i {
                tol = tol_i;
            }
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iters, false)
}

/// Multi-threading Bairstow's method (specific for auto-correlation function)
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
pub fn pbairstow_autocorr_mt(coeffs: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use threadpool::ThreadPool;

    let m_rs = vrs.len();
    let n_workers = 4; // assume 4 cores
    let (tx, rx) = channel();
    let pool = ThreadPool::new(n_workers);
    let mut converged = vec![false; m_rs];

    for niter in 1..options.max_iters {
        let mut tol = 0.0;
        let mut n_jobs = 0;

        for i in (0..m_rs).filter(|x| !converged[*x]) {
            let tx = tx.clone();
            let vrsc = vrs.clone();
            let mut pb = coeffs.to_owned();
            n_jobs += 1;
            pool.execute(move || {
                // let mut pb = coeffs.to_owned();
                let degree = pb.len() - 1; // assumed divided by 4
                let vri = vrsc[i];
                let mut vA = horner(&mut pb, degree, &vri);
                let tol_i = vA.norm_inf();
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                    return;
                }
                let mut vA1 = horner(&mut pb, degree - 2, &vri);
                for (_j, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
                    // vA1 -= delta(&vA, vrj, &(vri - vrj));
                    suppress_old(&mut vA, &mut vA1, &vri, vrj);
                    let vrjn = Vector2::<f64>::new(-vrj.x_, 1.0) / vrj.y_;
                    // vA1 -= delta(&vA, &vrjn, &(vri - vrjn));
                    suppress_old(&mut vA, &mut vA1, &vri, &vrjn);
                }
                let vrin = Vector2::<f64>::new(-vri.x_, 1.0) / vri.y_;
                // vA1 -= delta(&vA, &vrin, &(vri - vrin));
                suppress_old(&mut vA, &mut vA1, &vri, &vrin);
                let dt = delta(&vA, &vri, &vA1); // Gauss-Seidel fashion
                tx.send((Some((tol_i, dt)), i))
                    .expect("channel will be there waiting for a pool");
            });
        }
        for (res, i) in rx.iter().take(n_jobs) {
            if let Some(result) = res {
                let (toli, dt) = result;
                if tol < toli {
                    tol = toli;
                }
                vrs[i] -= dt;
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

/// Extract the quadratic function where its roots are within a unit circle
///
/// x^2 - r*x - t or x^2 + (r/t) * x + (-1/t)
/// (x - a1)(x - a2) = x^2 - (a1 + a2) x + a1 * a2
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
