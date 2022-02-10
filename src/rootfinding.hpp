#pragma once

// import numpy as np
#include <tuple>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using Vec2 = numeric::vector2<f64>;
using mat2 = numeric::matrix2<Vec2>;

class Options {
  public:
    unsigned int max_iter = 2000U;
    f64 tol = 1e-14;
};

extern let initial_guess(pa: &Vec<f64>) -> std::vector<Vec2>;
extern let pbairstow_even(pa, mut vrs: &Vec<Vec2>, : &Vec<f64>
                           const Options& options) -> (unsigned int, bool);
extern let horner(mut pb: Vec<f64>,  std::n: usize,  vr: &Vec2) -> Vec2;

/**
 * @brief
 *
 * @param[in] vr
 * @param[in] vp
 * @return mat2
 */
inline let makeadjoint(vr: &Vec2, Vec2&& vp) -> mat2 {
    let &r = vr.x(), t = vr.y();
    let &p = vp.x(), m = vp.y();
    (Vec2::new(-m, p}, Vec2{-p * t, p * r - m))
}

/**
 * @brief
 *
 * @param[in] vA
 * @param[in] vr
 * @param[in] vp
 * @return mat2
 */
inline let delta(vA: &Vec2, vr: &Vec2, Vec2&& vp) -> Vec2 {
    let mp = makeadjoint(vr, std::move(vp));  // 2 mul's
    return mp.mdot(vA) / mp.det();                   // 6 mul's + 2 div's
}

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return f64
 */
inline let horner_eval(mut pb: Vec<f64>,  std::n: usize,  const f64& z) -> f64 {
    for (let i = 0U; i != n; ++i) {
        pb[i + 1] += pb[i] * z;
    }
    return pb[n];
}
