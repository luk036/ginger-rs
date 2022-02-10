#include <bairstow/ThreadPool.h>  // for ThreadPool
#include <stddef.h>               // for usize

// #include <__bit_reference>           // for __bit_reference
#include <bairstow/rootfinding.hpp>  // for Vec2, delta, Options, horner_eval
#include <cmath>                     // for abs, acos, cos, pow
#include <functional>                // for __base
#include <future>                    // for future
#include <thread>                    // for thread
#include <tuple>                     // for tuple
#include <type_traits>               // for move
#include <vector>                    // for vector, vector<>::reference, __v...

#include "bairstow/vector2.hpp"  // for operator-, vector2

// using Vec2 = numeric::vector2<f64>;
// using mat2 = numeric::matrix2<Vec2>;

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] vr
 * @return Vec2
 */
pub fn horner(mut pb: Vec<f64>,  n: usize,  vr: &Vec2) -> Vec2 {
    let r = vr.x();
    let t = vr.y();
    pb[1] -= pb[0] * r;
    for (let i = 2U; i != n; ++i) {
        pb[i] -= pb[i - 1] * r + pb[i - 2] * t;
    }
    pb[n] -= pb[n - 2] * t;
    return Vec2::new(pb[n - 1], pb[n]);
}

/**
 * @brief
 *
 * @param[in] pa
 * @return Vec<Vec2>
 */
pub fn initial_guess(pa: &Vec<f64>) -> Vec<Vec2> {
    static let PI = std::acos(-1.0);

    let nn = pa.len() - 1;
    let c = -pa[1] / (nn * pa[0]);
    let pb = pa;
    let Pc = horner_eval(pb, nn, c);  // ???
    let re = std::pow(Pc.abs(), 1.0 / nn);
    nn /= 2;
    nn *= 2;  // make even
    let k = PI / nn;
    let m = c * c + re * re;
    let vr0s = Vec<Vec2>::new();
    for (let i = 1; i < nn; i += 2) {
        let temp = re * (k * i).cos();
        let r0 = -2 * (c + temp);
        let t0 = m + 2 * c * temp;
        vr0s.emplace_back(Vec2{std::move(r0), std::move(t0)});
    }
    return vr0s;
}

/**
 * @brief Multi-threading Bairstow's method (even degree only)
 *
 * @param[in] pa polynomial
 * @param[in,out] vrs vector of iterates
 * @param[in] options maximum iterations and tolorance
 * @return (unsigned int, bool)
 */
pub fn pbairstow_even(pa: &Vec<f64>, mut vrs: &Vec<Vec2>, 
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
                // let n = pa.len() - 1;
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
                }
                vrs[i] -= delta(vA, vri, std::move(vA1));  // Gauss-Seidel fashion
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
