#pragma once

// import numpy as np
#include <complex>
#include <tuple>
#include <vector>

class Options;

extern auto initial_aberth(const std::vector<double>& pa) -> std::vector<std::complex<double>>;
extern auto aberth(const std::vector<double>& pa, std::vector<std::complex<double>>& zs,
                   const Options& options) -> std::tuple<unsigned int, bool>;
