#pragma once

// import numpy as np
#include <tuple>
#include <vector>

#include "matrix2.hpp"
#include "vector2.hpp"

using Vec2 = numeric::vector2<f64>;
using mat2 = numeric::matrix2<Vec2>;

class Options;

extern let initial_autocorr(pa: &Vec<f64>) -> std::vector<Vec2>;
extern let pbairstow_autocorr(pa, mut vrs: &Vec<Vec2>, : &Vec<f64>
                               const Options& options) -> (unsigned int, bool);
extern void extract_autocorr(Vec2& vr);
