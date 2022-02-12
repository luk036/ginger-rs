#pragma once

// import numpy as np
#include <tuple>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using vec2 = numeric::vector2<double>;
using mat2 = numeric::matrix2<vec2>;

class Options {
  public:
    unsigned int max_iter = 2000U;
    double tol = 1e-14;
};

extern auto initial_guess(const std::vector<double>& pa) -> std::vector<vec2>;
extern auto pbairstow_even(const std::vector<double>& pa, std::vector<vec2>& vrs,
                           const Options& options) -> std::tuple<unsigned int, bool>;
extern auto horner(std::vector<double>& pb, std::size_t n, const vec2& vr) -> vec2;

/**
 * @brief
 *
 * @param[in] vr
 * @param[in] vp
 * @return mat2
 */
inline auto makeadjoint(const vec2& vr, vec2&& vp) -> mat2 {
    const auto &r = vr.x(), t = vr.y();
    const auto &p = vp.x(), m = vp.y();
    return {vec2{-m, p}, vec2{-p * t, p * r - m}};
}

/**
 * @brief
 *
 * @param[in] vA
 * @param[in] vr
 * @param[in] vp
 * @return mat2
 */
inline auto delta(const vec2& vA, const vec2& vr, vec2&& vp) -> vec2 {
    const auto mp = makeadjoint(vr, vp);  // 2 mul's
    return mp.mdot(vA) / mp.det();                   // 6 mul's + 2 div's
}

/**
 * @brief
 *
 * @param[in,out] pb
 * @param[in] n
 * @param[in] r
 * @return double
 */
inline auto horner_eval(std::vector<double>& pb, std::size_t n, const double& z) -> double {
    for (auto i = 0U; i != n; ++i) {
        pb[i + 1] += pb[i] * z;
    }
    return pb[n];
}
