// using Vec2 = numeric::vector2<f64>;
// using Mat2 = numeric::matrix2<Vec2>;

use std::sync::mpsc::channel;
use threadpool::ThreadPool;
// use rusty_pool::ThreadPool;
// use std::thread;
// use std::time::Duration;

use super::{Matrix2, Vector2};
type Vec2 = Vector2<f64>;
type Mat2 = Matrix2<f64>;

const PI: f64 = std::f64::consts::PI;

pub struct Options {
    pub max_iter: usize,
    pub tol: f64,
}

/**
 * @brief
 *
 * @param[in] vr
 * @param[in] vp
 * @return Mat2
 */
#[inline]
pub fn makeadjoint(vr: &Vec2, vp: &Vec2) -> Mat2 {
    let Vec2 { x_: r, y_: t } = *vr;
    let Vec2 { x_: p, y_: m } = *vp;
    Mat2::new(Vec2::new(-m, p), Vec2::new(-p * t, p * r - m))
}

/**
 * @brief
 *
 * @param[in] vaa
 * @param[in] vr
 * @param[in] vp
 * @return Mat2
 */
#[inline]
pub fn delta(vaa: &Vec2, vr: &Vec2, vp: &Vec2) -> Vec2 {
    let mp = makeadjoint(vr, vp); // 2 mul's
    mp.mdot(vaa) / mp.det() // 6 mul's + 2 div's
}

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return f64
 */
#[inline]
pub fn horner_eval(pb: &mut [f64], n: usize, z: f64) -> f64 {
    for i in 0..n {
        pb[i + 1] += pb[i] * z;
    }
    pb[n]
}

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] vr
 * @return Vec2
 */
pub fn horner(pb: &mut [f64], n: usize, vr: &Vec2) -> Vec2 {
    let Vec2 { x_: r, y_: t } = vr;
    pb[1] -= pb[0] * r;
    for i in 2..n {
        pb[i] -= pb[i - 1] * r + pb[i - 2] * t;
    }
    pb[n] -= pb[n - 2] * t;
    Vec2::new(pb[n - 1], pb[n])
}

/**
 * @brief
 *
 * @param[in] pa
 * @return Vec<Vec2>
 */
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
        vr0s.push(Vec2::new(r0, t0));
    }
    vr0s
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param[in] pa polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (usize, bool)
 */
pub fn pbairstow_even(pa: &[f64], vrs: &mut Vec<Vec2>, options: &Options) -> (usize, bool) {
    // let n = pa.len() - 1; // degree, assume even
    let m = vrs.len();
    let mut found = false;
    let mut converged = vec![false; m];

    let mut niter: usize = 0;
    while niter != options.max_iter {
        niter += 1;

        let mut tol = 0.0;
        let mut rx = vec![];

        for i in 0..m {
            if converged[i] {
                continue;
            }
            let mut job = || {
                let mut pb = pa.to_owned();
                let n = pa.len() - 1; // degree, assume even
                let m = vrs.len();
                let vri = vrs[i];
                let vaa = horner(&mut pb, n, &vri);
                let tol_i = vaa.norm_inf();
                if tol_i < 1e-15 {
                    converged[i] = true;
                    rx.push(tol_i);
                } else {
                    let mut vaa1 = horner(&mut pb, n - 2, &vri);
                    for j in 0..m {
                        // exclude i
                        if j == i {
                            continue;
                        }
                        let vrj = vrs[j].to_owned(); // make a copy, don't reference!
                        vaa1 -= delta(&vaa, &vrj, &(vri - vrj));
                    }
                    vrs[i] -= delta(&vaa, &vri, &vaa1); // Gauss-Seidel fashion
                                                        // tx.send(tol_i).expect("channel will be there waiting for a pool");
                    rx.push(tol_i);
                    // return tol_i;
                }
            };

            job();
            // }));
        }
        for result in rx.iter() {
            if tol < *result {
                tol = *result;
            }
        }
        if tol < options.tol {
            found = true;
            break;
        }
    }
    (niter, found)
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param[in] pa polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (usize, bool)
 */
pub fn pbairstow_even_th(pa: Vec<f64>, vrs: Vec<Vec2>, options: &Options) -> (usize, bool) {
    let m = vrs.len();
    let mut found = false;
    let n_workers = 4; // assume 4 cores
    // let mut converged = vec![0_u32; m];

    let mut niter: usize = 0;
    while niter != options.max_iter {
        niter += 1;

        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);

        for i in 0..m {
            // if converged[i] {
            //     continue;
            // }
            let tx = tx.clone();
            let mut vrsc = vrs.clone();
            let mut pb = pa.clone();
            pool.execute(move || {
                // let mut pb = pa.to_owned();
                let n = pb.len() - 1; // degree, assume even
                let m = vrsc.len();
                let vri = vrsc[i];
                let vaa = horner(&mut pb, n, &vri);
                let tol_i = vaa.norm_inf();
                if tol_i < 1e-15 {
                    // *&mut converged[i] += 1;
                    tx.send(None)
                        .expect("channel will be there waiting for a pool");
                } else {
                    let mut vaa1 = horner(&mut pb, n - 2, &vri);
                    for j in 0..m {
                        // exclude i
                        if j == i {
                            continue;
                        }
                        let vrj = vrsc[j]; // make a copy, don't reference!
                        vaa1 -= delta(&vaa, &vrj, &(vri - vrj));
                    }
                    vrsc[i] -= delta(&vaa, &vri, &vaa1); // Gauss-Seidel fashion
                    tx.send(Some(tol_i))
                        .expect("channel will be there waiting for a pool");
                }
            });
        }
        for res in rx.iter() {
            // let result = (*job).await_complete();
            if let Some(result) = res {
                if tol < result {
                    tol = result;
                }
            } else {
                // don't care
            } 
        }
        if tol < options.tol {
            found = true;
            break;
        }
    }
    (niter, found)
}

// let find_rootq(r: &Vec2) {
//     let hb = b / 2.;
//     let d = hb * hb - c;
//     if d < 0.  {
//         let x1 = -hb + (sqrt(-d) if hb < 0. else -sqrt(-d )*1j;
//     }
//     else {
//         let x1 = -hb + (d.sqrt() if hb < 0. else -sqrt(d );
//     }
//     let x2 = c / x1;
//     return x1, x2;
// }
