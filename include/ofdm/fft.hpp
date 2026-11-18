#pragma once

#include <cstddef>

#include "ofdm/types.hpp"

namespace ofdm {

bool isPowerOfTwo(std::size_t n);

// Iterative radix-2 Cooley-Tukey FFT, O(N log N). Size must be a power of two.
// Conventions (same as MATLAB/NumPy):
//   fft:  X[k] = sum_n x[n] e^{-j2pi kn/N}
//   ifft: x[n] = (1/N) sum_k X[k] e^{+j2pi kn/N}
CVec fft(CVec x);
CVec ifft(CVec x);

// Direct O(N^2) transforms, any size. Kept as a reference to validate the FFT.
CVec dft(const CVec& x);
CVec idft(const CVec& x);

}  // namespace ofdm
