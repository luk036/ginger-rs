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
    let ans = Complex::<f64>::new(pb[0], 0.0);
    for i in 1..pb.len() {
        ans = ans * z + pb[i];
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
pub fn horner_eval_f(pb: &[f64], z: f64) -> f64 {
    let ans = pb[0]
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
    let Pc = horner_eval_f(pa, c);
    let re = (-Pc as Complex<f64>).pow(1.0 / n); // ???
    let k = TWO_PI / (n as f64);
    let z0s = vec![];
    for i in 0..n {
        let theta = k * (i + 0.25);
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
 * @return (unsigned int, bool)
 */
pub fn aberth(pa, std::vector<std::complex<f64>>& zs,: &Vec<f64>
            options: &Options) -> (unsigned int, bool) {
    let m = zs.len();
    let n = pa.len() - 1;  // degree, assume even
    let found = false;
    let converged = vec![false; m];
    let pb = vec![0.0; n];
    for i in 0..n {
        pb[i] = (n - i) * pa[i];
    }
    // let niter = 1U;
    // ThreadPool pool(std::thread::hardware_concurrency());

    for niter in 1..options.max_iter {
        let tol = 0.0;
        // std::vector<std::future<f64>> results;

        for i in 0..m {
            if converged[i]  {
                continue;
            }
            // results.emplace_back(pool.enqueue([&, i]() {
                let zi = zs[i];
                let P = horner_eval_g(pa, zi);
                let tol_i = P.abs();
                if tol_i < 1e-15  {
                    converged[i] = true;
                    return tol_i;
                }
                let P1 = horner_eval_g(pb, zi);
                for j in 0..m {  // exclude i
                    if j == i  {
                        continue;
                    }
                    let zj = zs[j];  // make a copy, don't reference!
                    P1 -= P / (zi - zj);
                }
                zs[i] -= P / P1;  // Gauss-Seidel fashion
                return tol_i;
            }));
        }
        for result : results  {
            let& res = result.get();
            if tol < res  {
                tol = res;
            }
        }
        if tol < options.tol  {
            found = true;
            break;
        }
    }
    (niter, found)
}
