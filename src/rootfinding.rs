use super::{Matrix2, Vector2};

type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;

const PI: f64 = std::f64::consts::PI;

pub struct Options {
    pub max_iter: usize,
    pub tol: f64,
    pub tol_ind: f64,
}

impl Default for Options {
    fn default() -> Self {
        Options {
            max_iter: 2000,
            tol: 1e-12,
            tol_ind: 1e-15,
        }
    }
}

/**
 * @brief
 *
 * @param vr
 * @param vp
 * @return Mat2
 */
#[inline]
pub fn makeadjoint(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let (r, t) = (vr.x_, vr.y_);
    let (p, m) = (vp.x_, vp.y_);
    Mat2::new(
        Vector2::<f64>::new(-m, p),
        Vector2::<f64>::new(-p * t, p * r - m),
    )
}

/// Extract the quadratic function where its roots are within a unit circle
///
/// r * p - m   -p
/// q * p       -m
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::delta;
/// use bairstow::vector2::Vector2;
///
/// let vd = delta(&Vector2::new(1.0, 2.0), &Vector2::new(2.0, 0.0), &Vector2::new(4.0, 5.0));
///
/// assert_eq!(vd, Vector2::new(-0.2, -0.4));
/// ```
#[inline]
pub fn delta(vaa: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let mp = makeadjoint(vr, vp); // 2 mul's
    mp.mdot(vaa) / mp.det() // 6 mul's + 2 div's
}

/// For ri - rj
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::delta1;
/// use bairstow::vector2::Vector2;
///
/// let vd = delta1(&Vector2::new(1.0, 2.0), &Vector2::new(-2.0, 0.0), &Vector2::new(4.0, -5.0));
///
/// assert_eq!(vd, Vector2::new(0.2, 0.4));
/// ```
#[inline]
pub fn delta1(vaa: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let (r, q) = (vr.x_, vr.y_);
    let (p, s) = (vp.x_, vp.y_);
    let mp = Matrix2::new(Vec2::new(-s, -p), Vec2::new(p * q, p * r - s));
    mp.mdot(vaa) / mp.det() // 6 mul's + 2 div's
}

#[inline]
pub fn suppress_old(vA: &mut Vec2, vA1: &mut Vec2, vri: &Vec2, vrj: &Vec2) {
    let (A, B) = (vA.x_, vA.y_);
    let (A1, B1) = (vA1.x_, vA1.y_);
    let vp = vri - vrj;
    let (r, q) = (vri.x_, -vri.y_);
    let (p, s) = (vp.x_, -vp.y_);
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

/// Horner evalution
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::horner_eval;
/// use approx_eq::assert_approx_eq;
///
/// let mut pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner_eval(&mut pa, 8, 2.0);
///
/// assert_approx_eq!(px, 18250.0);
/// assert_approx_eq!(pa[3], 460.0);
/// ```
#[inline]
pub fn horner_eval(pb: &mut [f64], n: usize, z: f64) -> f64 {
    for i in 0..n {
        pb[i + 1] += pb[i] * z;
    }
    pb[n]
}

/// Horner evalution for Bairstow's method
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::horner;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let mut pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let px = horner(&mut pa, 8, &Vector2::new(1.0, 2.0));
///
/// assert_approx_eq!(px.x_, 114.0);
/// assert_approx_eq!(px.y_, 134.0);
/// assert_approx_eq!(pa[3], 15.0);           
/// ```
pub fn horner(pb: &mut [f64], n: usize, vr: &Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: t } = vr;
    pb[1] -= pb[0] * r;
    for i in 2..n {
        pb[i] -= pb[i - 1] * r + pb[i - 2] * t;
    }
    pb[n] -= pb[n - 2] * t;
    Vector2::<f64>::new(pb[n - 1], pb[n])
}

/// Initial guess for Bairstow's method
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::initial_guess;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_guess(&pa);
///
/// assert_approx_eq!(vr0s[0].x_, -1.4537520283010816);
/// assert_approx_eq!(vr0s[0].y_, 0.7559947795353074);
/// ```
pub fn initial_guess(pa: &[f64]) -> Vec<Vec2> {
    let mut n = pa.len() - 1;
    let c = -pa[1] / (pa[0] * n as f64);
    let mut pb = pa.to_owned();
    let centroid = horner_eval(&mut pb, n, c); // ???
    let re = centroid.abs().powf(1.0 / (n as f64));
    n /= 2;
    n *= 2; // make even
    let k = PI / (n as f64);
    let m = c * c + re * re;
    let mut vr0s = Vec::<Vec2>::new();
    for i in (1..n).step_by(2) {
        let temp = re * (k * i as f64).cos();
        let r0 = -2.0 * (c + temp);
        let t0 = m + 2.0 * c * temp;
        vr0s.push(Vector2::<f64>::new(r0, t0));
    }
    vr0s
}

/// Bairstow's method (even degree only)
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_guess, pbairstow_even, Options};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&pa);
/// let (niter, _found) = pbairstow_even(&pa, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 9);
/// ```
pub fn pbairstow_even(pa: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    let n = pa.len() - 1; // degree, assume even
    let m = vrs.len();
    let mut converged = vec![false; m];

    for niter in 1..options.max_iter {
        let mut tol = 0.0;
        for i in 0..m {
            if converged[i] {
                continue;
            }
            let mut pb = pa.to_owned();
            let vri = vrs[i];
            let vaa = &horner(&mut pb, n, &vri);
            let tol_i = vaa.norm_inf();
            if tol_i < 1e-15 {
                converged[i] = true;
                continue;
            } else {
                let mut vaa1 = horner(&mut pb, n - 2, &vri);
                if tol < tol_i {
                    tol = tol_i;
                }
                for (j, vrj) in vrs.iter().enumerate().take(m) {
                    if j == i {
                        continue;
                    }
                    vaa1 -= delta(vaa, vrj, &(vri - vrj));
                }
                vrs[i] -= delta(vaa, &vri, &vaa1); // Gauss-Seidel fashion
            }
        }
        if tol < options.tol {
            return (niter, true);
        }
    }
    (options.max_iter, false)
}

/// Multi-threading Bairstow's method (even degree only)
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_guess, pbairstow_even_th, Options};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_guess(&pa);
/// let (niter, _found) = pbairstow_even_th(&pa, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 7);
/// ```
pub fn pbairstow_even_th(pa: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use threadpool::ThreadPool;

    let n_workers = 4; // assume 4 cores
    let m = vrs.len();
    let mut converged = vec![false; m];

    for niter in 1..options.max_iter {
        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);
        let mut n_jobs = 0;
        for i in (0..m).filter(|x| !converged[*x]) {
            // if converged[i] {
            //     continue;
            // }
            let tx = tx.clone();
            let vrsc = vrs.clone();
            let mut pb = pa.to_owned();

            n_jobs += 1;
            pool.execute(move || {
                // let mut pb = pa.to_owned();
                let n = pb.len() - 1; // degree, assume even
                let vri = vrsc[i];
                let vaa = horner(&mut pb, n, &vri);
                let tol_i = vaa.norm_inf();
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                } else {
                    let mut vaa1 = horner(&mut pb, n - 2, &vri);
                    for (_, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
                        vaa1 -= delta(&vaa, vrj, &(vri - vrj));
                    }
                    let dt = delta(&vaa, &vri, &vaa1); // Gauss-Seidel fashion
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
    (options.max_iter, false)
}

/// Initial guess for Bairstow's method (specific for auto-correlation function)
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::initial_autocorr;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let vr0s = initial_autocorr(&pa);
///
/// assert_approx_eq!(vr0s[0].x_, -1.8858840950805662);
/// assert_approx_eq!(vr0s[0].y_, 1.7782794100389228);
/// ```
pub fn initial_autocorr(pa: &[f64]) -> Vec<Vec2> {
    let mut n = pa.len() - 1;
    let re = (pa[n].abs() as f64).powf(1.0 / (n as f64));
    n /= 2;
    let k = PI / (n as f64);
    let m = re * re;
    let mut vr0s = Vec::<Vec2>::new();
    for i in (1..n).step_by(2) {
        vr0s.push(Vector2::<f64>::new(-2.0 * re * (k * i as f64).cos(), m));
    }
    vr0s
}

/// Simultenous Bairstow's method (specific for auto-correlation function)
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_autocorr, pbairstow_autocorr, Options};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&pa);
/// let (niter, _found) = pbairstow_autocorr(&pa, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn pbairstow_autocorr(pa: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    let m = vrs.len();
    let mut converged = vec![false; m];

    for niter in 1..options.max_iter {
        let mut tol = 0.0;

        for i in 0..m {
            if converged[i] {
                continue;
            }
            let mut job = || {
                let mut pb = pa.to_owned();
                let n = pa.len() - 1; // assumed divided by 4
                let vri = vrs[i];
                let vaa = horner(&mut pb, n, &vri);
                let tol_i = vaa.norm_inf();
                if tol_i < 1e-15 {
                    converged[i] = true;
                    return tol_i;
                }
                let mut vaa1 = horner(&mut pb, n - 2, &vri);
                for (_j, vrj) in vrs.iter().enumerate().filter(|t| t.0 != i) {
                    vaa1 -= delta(&vaa, vrj, &(vri - vrj));
                    let vrjn = Vector2::<f64>::new(vrj.x_, 1.0) / vrj.y_;
                    vaa1 -= delta(&vaa, &vrjn, &(vri - vrjn));
                }
                let vrin = Vector2::<f64>::new(vri.x_, 1.0) / vri.y_;
                vaa1 -= delta(&vaa, &vrin, &(vri - vrin));
                vrs[i] -= delta(&vaa, &vri, &vaa1); // Gauss-Seidel fashion
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
    (options.max_iter, false)
}

/// Multi-threading Bairstow's method (specific for auto-correlation function)
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::{initial_autocorr, pbairstow_autocorr_th, Options};
///
/// let pa = vec![10.0, 34.0, 75.0, 94.0, 150.0, 94.0, 75.0, 34.0, 10.0];
/// let mut vrs = initial_autocorr(&pa);
/// let (niter, _found) = pbairstow_autocorr_th(&pa, &mut vrs, &Options::default());
///
/// assert_eq!(niter, 5);
/// ```
pub fn pbairstow_autocorr_th(pa: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use threadpool::ThreadPool;

    let m = vrs.len();
    let n_workers = 4; // assume 4 cores
    let (tx, rx) = channel();
    let pool = ThreadPool::new(n_workers);
    let mut converged = vec![false; m];

    for niter in 1..options.max_iter {
        let mut tol = 0.0;
        let mut n_jobs = 0;

        for i in (0..m).filter(|x| !converged[*x]) {
            let tx = tx.clone();
            let vrsc = vrs.clone();
            let mut pb = pa.to_owned();
            n_jobs += 1;
            pool.execute(move || {
                // let mut pb = pa.to_owned();
                let n = pb.len() - 1; // assumed divided by 4
                let vri = vrsc[i];
                let vaa = horner(&mut pb, n, &vri);
                let tol_i = vaa.norm_inf();
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                    return;
                }
                let mut vaa1 = horner(&mut pb, n - 2, &vri);
                for (_j, vrj) in vrsc.iter().enumerate().filter(|t| t.0 != i) {
                    vaa1 -= delta(&vaa, vrj, &(vri - vrj));
                    let vrjn = Vector2::<f64>::new(vrj.x_, 1.0) / vrj.y_;
                    vaa1 -= delta(&vaa, &vrjn, &(vri - vrjn));
                }
                let vrin = Vector2::<f64>::new(vri.x_, 1.0) / vri.y_;
                vaa1 -= delta(&vaa, &vrin, &(vri - vrin));
                let dt = delta(&vaa, &vri, &vaa1); // Gauss-Seidel fashion
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
    (options.max_iter, false)
}

/// Extract the quadratic function where its roots are within a unit circle
///
/// x^2 + r*x + t or x^2 + (r/t) * x + (1/t)
/// (x + a1)(x + a2) = x^2 + (a1 + a2) x + a1 * a2
///
/// Examples:
///
/// ```
/// use bairstow::rootfinding::extract_autocorr;
/// use bairstow::vector2::Vector2;
/// use approx_eq::assert_approx_eq;
///
/// let vr = extract_autocorr(Vector2::new(-1.0, 4.0));
///
/// assert_approx_eq!(vr.x_, -0.25);
/// assert_approx_eq!(vr.y_, 0.25);
/// ```
#[allow(dead_code)]
pub fn extract_autocorr(vr: Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: t } = vr;
    let hr = r / 2.0;
    let d = hr * hr - t;
    if d < 0.0 {
        // complex conjugate root
        if t > 1.0 {
            return Vector2::<f64>::new(r, 1.0) / t;
        }
    }
    // two real roots
    let mut a1 = hr + (if hr >= 0.0 { d.sqrt() } else { -d.sqrt() });
    let mut a2 = t / a1;

    if a1.abs() > 1.0 {
        if a2.abs() > 1.0 {
            a2 = 1.0 / a2;
        }
        a1 = 1.0 / a1;
        return Vector2::<f64>::new(a1 + a2, a1 * a2);
    }
    if a2.abs() > 1.0 {
        a2 = 1.0 / a2;
        return Vector2::<f64>::new(a1 + a2, a1 * a2);
    }
    // else no need to change
    vr
}
