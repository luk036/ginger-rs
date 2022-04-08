use super::Options;
use num::Complex;
// use lds_rs::lds::Circle;

const TWO_PI: f64 = std::f64::consts::TAU;

/**
 * @brief
 *
 * @param pb
 * @param n
 * @param r
 * @return f64
 */
pub fn horner_eval_c(pb: &[f64], z: &Complex<f64>) -> Complex<f64> {
    let mut ans = Complex::<f64>::new(pb[0], 0.0);
    for coeff in pb.iter().skip(1) {
        ans *= z;
        ans += coeff;
    }
    ans
}

/**
 * @brief
 *
 * @param pb
 * @param n
 * @param r
 * @return f64
 */
pub fn horner_eval_f(pb: &[f64], z: &f64) -> f64 {
    let mut ans = pb[0];
    for coeff in pb.iter().skip(1) {
        ans = ans * z + coeff;
    }
    ans
}

/**
 * @brief
 *
 * @param pa
 * @return Vec<Vec2>
 */
pub fn initial_aberth(pa: &[f64]) -> Vec<Complex<f64>> {
    let n = pa.len() - 1;
    let c = -pa[1] / (pa[0] * n as f64);
    let ppc = horner_eval_f(pa, &c);
    let re = Complex::<f64>::new(-ppc, 0.0).powf(1.0 / n as f64);
    let k = TWO_PI / (n as f64);
    let mut z0s = vec![];
    for i in 0..n {
        let theta = k * (0.25 + i as f64);
        let z0 = c + re * Complex::<f64>::new(theta.cos(), theta.sin());
        z0s.push(z0);
    }
    z0s
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param pa polynomial
 * @param zs vector of iterates
 * @param options maximum iterations and tolorance
 * @return (usize, bool)
 */
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

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param pa polynomial
 * @param zs vector of iterates
 * @param options maximum iterations and tolorance
 * @return (usize, bool)
 */
pub fn aberth_th(pa: &[f64], zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    use std::sync::mpsc::channel;
    use std::sync::Arc;
    use threadpool::ThreadPool;

    let n_workers = 4; // assume 4 cores

    let m = zs.len();
    let n = pa.len() - 1; // degree, assume even
    let mut pb = vec![0.0; n];
    let n = pa.len() - 1; // degree, assume even
    for k in 0..n {
        pb[k] = pa[k] * (n - k) as f64;
    }
    // let mut zsc = zs.clone();
    let pb = pb; // make imutatable
    let pa_share = Arc::new(pa.to_owned());
    let pb_share = Arc::new(pb);
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
