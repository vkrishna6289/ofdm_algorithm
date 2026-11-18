#pragma once

#include <complex>
#include <vector>

namespace ofdm {

using Complex = std::complex<double>;
using CVec = std::vector<Complex>;
using Bits = std::vector<int>;  // each element is 0 or 1

constexpr double kPi = 3.14159265358979323846;

// |z|^2. Use this instead of std::norm in hot loops: libstdc++ implements
// std::norm as abs(z)^2, i.e. via hypot(), which is very slow on MinGW.
inline double mag2(const Complex& z) { return z.real() * z.real() + z.imag() * z.imag(); }

}  // namespace ofdm
