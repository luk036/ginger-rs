const PI: f64 = std::f64::consts::PI;

/**
 * @brief initial guess (specific for let-correlation function)
 *
 * @param[in] pa
 * @return Vec<Vec2>
 */
pub fn initial_autocorr(pa: &Vec<f64>) -> Vec<Vec2> {
    let nn = pa.len() - 1;
    let re = (pa[nn].abs() as f64).powf(1.0 / (nn as f64));

    nn /= 2;
    let k = PI / (nn as f64);
    let m = re * re;
    let vr0s = Vec<Vec2>::new();
    for i in (1..nn).step_by(2) {
        vr0s.push(Vec2::new(-2.0 * re * (k * i as f64).cos(), m));
    }
    return vr0s;
}

/**
 * @brief Multi-threading Bairstow's method (specific for let-correlation function)
 *
 * @param[in] pa polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (unsigned int, bool)
 */
pub fn pbairstow_autocorr(pa: &Vec<f64>, mut vrs: &Vec<Vec2>, 
                        options: &Options) -> (unsigned int, bool) {
    let nn = pa.len() - 1;  // degree, assume even
    let mm = vrs.len();
    let found = false;
    let converged = vec![false; mm];
    // let niter = 1U;
    // ThreadPool pool(std::thread::hardware_concurrency());

    for niter in 1..options.max_iter {
        let tol = 0.0;
        // std::vector<std::future<f64>> results;
        for i in 0..mm {
            if converged[i]  {
                continue;
            }
            // results.emplace_back(pool.enqueue([&, i]() {
                let pb = pa;
                let vri = vrs[i];
                let vA = horner(pb, nn, vri);
                let tol_i = vA.norm_inf();
                if tol_i < 1e-15  {
                    converged[i] = true;
                    return tol_i;
                }
                let vA1 = horner(pb, nn - 2, vri);
                for j in 0..mm {  // exclude i
                    if j == i  {
                        continue;
                    }
                    let vrj = vrs[j];  // make a copy, don't reference!
                    vA1 -= delta(vA, vrj, vri - vrj);
                    let vrjn = numeric::vector2<f64>(vrj.x(), 1.0) / vrj.y();
                    vA1 -= delta(vA, vrjn, vri - vrjn);
                }
                let vrin = numeric::vector2<f64>(vri.x(), 1.0) / vri.y();
                vA1 -= delta(vA, vrin, vri - vrin);

                vrs[i] -= delta(vA, vri, vA1);  // Gauss-Seidel fashion
                return tol_i;
            }));
        }
        for result : results  {
            let& res = result.get();
            if tol < res  {
                tol = res;
            }
        }
        // fmt::print("tol: {}\n", tol);
        if tol < options.tol  {
            found = true;
            break;
        }
    }
    (niter, found)
}

/**
 * @brief Extract the quadratic function where its roots are within a unit circle
 *
 *   x^2 + r*x + t or x^2 + (r/t) * x + (1/t)
 *   (x + a1)(x + a2) = x^2 + (a1 + a2) x + a1 * a2
 *
 * @param[in,out] vr
 */
void extract_autocorr(Vec2& vr) {
    let r = vr.x();
    let t = vr.y();
    let hr = r / 2.0;
    let d = hr * hr - t;
    if d < 0.0  {  // complex conjugate root
        if t > 1.0  {
            vr = Vec2::new(r, 1.0) / t;
        }
        // else no need to change
    } else {  // two real roots
        let a1 = hr + (if hr >= 0.0 { d.sqrt() } else { -d.sqrt() });
        let a2 = t / a1;
        if a1.abs() > 1.0 {
            if a2.abs() > 1.0 {
                a2 = 1.0 / a2;
            }
            a1 = 1.0 / a1;
            vr = Vec2::new(a1 + a2, a1 * a2);
        } else if a2.abs() > 1.0 {
            a2 = 1.0 / a2;
            vr = Vec2::new(a1 + a2, a1 * a2);
        }
        // else no need to change
    }
}
