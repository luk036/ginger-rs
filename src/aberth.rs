use super::Options;
use num::Complex;

const TWO_PI: f64 = std::f64::consts::TAU;

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return f64
 */
pub fn horner_eval_c(pb: &[f64], z: &Complex<f64>) -> Complex<f64> {
    let mut ans = Complex::<f64>::new(pb[0], 0.0);
    for i in 1..pb.len() {
        ans *= z;
        ans += pb[i];
    }
    ans
}

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return f64
 */
pub fn horner_eval_f(pb: &[f64], z: &f64) -> f64 {
    let mut ans = pb[0];
    for i in 1..pb.len() {
        ans = ans * z + pb[i];
    }
    ans
}

/**
 * @brief
 *
 * @param[in] pa
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
 * @param[in] pa polynomial
 * @param[in,out] zs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (usize, bool)
 */
pub fn aberth(pa: &Vec<f64>, zs: &mut Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
    let m = zs.len();
    let n = pa.len() - 1; // degree, assume even
    let mut found = false;
    let mut converged = vec![false; m];
    let mut pb = vec![0.0; n];
    for i in 0..n {
        pb[i] = pa[i] * (n - i) as f64;
    }

    let mut niter: usize = 0;
    while niter < options.max_iter {
        niter += 1;

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
                for j in 0..m {
                    // exclude i
                    if j == i {
                        continue;
                    }
                    let zj = zs[j]; // make a copy, don't reference!
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
 * @param[in,out] zs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (usize, bool)
 */
pub fn aberth_th(pa: Vec<f64>, mut zs: Vec<Complex<f64>>, options: &Options) -> (usize, bool) {
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
    let pb = pb; // make imutatable
    let pa_share = Arc::new(pa);
    let pb_share = Arc::new(pb);
    let zs_share = Arc::new(zs);

    let mut found = false;
    let mut converged = vec![false; m];

    let mut niter: usize = 0;
    while niter < options.max_iter {
        niter += 1;

        let mut tol = 0.0;
        let (tx, rx) = channel();
        let pool = ThreadPool::new(n_workers);
        let mut n_jobs = 0;

        for i in 0..m {
            if converged[i] {
                continue;
            }
            let tx = tx.clone();
            // let zi = Complex::<f64>::default();
            let pa_clone = Arc::clone(&pa_share);
            let pb_clone = Arc::clone(&pb_share);
            let zs_clone = Arc::clone(&zs_share);
    
            // let zs_ref = &zs;
            // let pa_ref = &pa;
            // let pb_ref = &pb;
            n_jobs += 1;
            pool.execute(move || {
                let zi = zs_clone[i];
                let pp = horner_eval_c(&pa_clone, &zi);
                let tol_i = pp.l1_norm(); // ???
                if tol_i < 1e-15 {
                    tx.send((None, i))
                        .expect("channel will be there waiting for a pool");
                }
                let mut pp1 = horner_eval_c(&pb_clone, &zi);
                for (j, zj) in zs_clone.iter().enumerate() {
                    // exclude i
                    if j == i {
                        continue;
                    }
                    pp1 -= pp / (zi - zj);
                }
                // zs[i] -= pp / pp1; // Gauss-Seidel fashion
                tx.send((Some(tol_i), i))
                    .expect("channel will be there waiting for a pool");
            });
        }
        for (res, i) in rx.iter().take(n_jobs) {
            if let Some(result) = res {
                let toli = result;
                if tol < toli {
                    tol = toli;
                }
                // vrs[i] -= dt;
            } else {
                converged[i] = true;
            }
        }
        if tol < options.tol {
            found = true;
            break;
        }
    }
    (niter, found)
}
