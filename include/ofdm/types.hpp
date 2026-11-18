#pragma once

#include <complex>
#include <vector>

namespace ofdm {

using Complex = std::complex<double>;
using CVec = std::vector<Complex>;
using Bits = std::vector<int>;  // each element is 0 or 1

constexpr double kPi = 3.14159265358979323846;

}  // namespace ofdm
