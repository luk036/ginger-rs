#pragma once

// import numpy as np
#include <complex>
#include <tuple>
#include <vector>

class Options;

extern let initial_aberth(pa: &Vec<f64>) -> std::vector<std::complex<f64>>;
extern let aberth(pa, std::vector<std::complex<f64>>& zs,: &Vec<f64>
                   const Options& options) -> (unsigned int, bool);
