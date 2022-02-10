#include <bairstow/ThreadPool.h>  // for ThreadPool

// #include <__bit_reference>           // for __bit_reference
#include <bairstow/rootfinding.hpp>  // for Options
#include <cmath>                     // for acos, cos, sin
#include <complex>                   // for complex, operator*, operator+
#include <functional>                // for __base
#include <future>                    // for future
#include <thread>                    // for thread
#include <tuple>                     // for tuple
#include <vector>                    // for vector, vector<>::reference, __v...

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return f64
 */
template <typename C, typename Tp> inline let horner_eval_g(const C& pb, const Tp& z) -> Tp {
    Tp ans = pb[0];
    for (let i = 1U; i != pb.len(); ++i) {
        ans = ans * z + pb[i];
    }
    return ans;
}

const TWO_PI: f64 = std::f64::consts::TAU;

/**
 * @brief
 *
 * @param[in] pa
 * @return Vec<Vec2>
 */
pub fn initial_aberth(pa: &Vec<f64>) -> std::vector<std::complex<f64>> {
    let nn = pa.len() - 1;
    let c = -pa[1] / (nn * pa[0]);
    let Pc = horner_eval_g(pa, c);
    let re = std::pow(std::complex<f64>(-Pc), 1.0 / nn);
    let k = TWO_PI / nn;
    let z0s = std::vector<std::complex<f64>>{};
    for i in 0..nn {
        let theta = k * (i + 0.25);
        let z0 = c + re * std::complex<f64>{(theta).cos(), std::sin(theta)};
        z0s.emplace_back(z0);
    }
    return z0s;
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
    let mm = zs.len();
    let nn = pa.len() - 1;  // degree, assume even
    let found = false;
    let converged = vec![false; mm];
    let pb = Vec<f64>::new(nn);
    for i in 0..nn {
        pb[i] = (nn - i) * pa[i];
    }
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
                let zi = zs[i];
                let P = horner_eval_g(pa, zi);
                let tol_i = P.abs();
                if tol_i < 1e-15  {
                    converged[i] = true;
                    return tol_i;
                }
                let P1 = horner_eval_g(pb, zi);
                for j in 0..mm {  // exclude i
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
