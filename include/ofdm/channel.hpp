#pragma once

#include <cstddef>
#include <random>

#include "ofdm/types.hpp"

namespace ofdm {

// Adds circularly-symmetric complex Gaussian noise with total variance
// noiseVariance (noiseVariance/2 per real dimension).
CVec addAwgn(const CVec& x, double noiseVariance, std::mt19937_64& rng);

// Linear convolution with a FIR channel, y[n] = sum_l h[l] x[n-l].
// Output has the same length as x (the tail past the last sample is dropped).
CVec applyMultipath(const CVec& x, const CVec& taps);

// Random Rayleigh-fading tap set with an exponential power-delay profile:
// E|h[l]|^2 proportional to exp(-l / decay), normalised so sum_l E|h[l]|^2 = 1.
// The channel is then frequency-selective but has unit average gain.
CVec rayleighTaps(std::size_t numTaps, double decay, std::mt19937_64& rng);

// Channel frequency response on nFft bins (FFT of the zero-padded taps).
CVec frequencyResponse(const CVec& taps, std::size_t nFft);

}  // namespace ofdm
