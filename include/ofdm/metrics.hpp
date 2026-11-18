#pragma once

#include <cmath>
#include <cstddef>
#include <vector>

#include "ofdm/ofdm.hpp"
#include "ofdm/types.hpp"

namespace ofdm {

// mean |x[n]|^2
double meanPower(const CVec& x);

// Welch power spectral density estimate: Hann-windowed segments of
// segmentLength samples, 50% overlap, averaged |FFT|^2 normalised by the
// window energy (white noise of variance s2 gives a flat PSD of s2).
// Returned in FFT bin order (bin 0 = DC). segmentLength must be a power of two.
// With segmentLength == p.nFft, bin k corresponds to subcarrier k.
std::vector<double> welchPsd(const CVec& x, std::size_t segmentLength);

// Sum of PSD over logical bins kLow..kHigh inclusive (negative k wrap to the top).
double bandPower(const std::vector<double>& psd, int kLow, int kHigh);

// RMS error vector magnitude relative to the reference power:
// sqrt( sum|r - ref|^2 / sum|ref|^2 ).  Multiply by 100 for percent.
double evmRms(const CVec& received, const CVec& reference);

// PAPR = max|x[n]|^2 / mean|x[n]|^2, in dB.
double paprDb(const CVec& x);

// Time-domain OFDM symbol (no CP) oversampled by factor L, obtained by
// zero-padding the middle of the spectrum and taking an L*N-point IFFT.
// Nyquist-rate samples (L = 1) miss the true analogue peaks; L >= 4 is the
// usual choice for accurate PAPR measurement.
CVec oversampledSymbol(const CVec& bins, std::size_t L);

std::size_t countBitErrors(const Bits& a, const Bits& b);

// Noise variance per time-domain sample that gives the requested Eb/N0 per
// data subcarrier for QPSK, given that the transmitter uses ifft() (1/N scaling)
// and data symbols have unit energy:
//   after the FFT, noise per bin = N * sigma^2 and Es = 1, Es/N0 = 2 Eb/N0
//   => sigma^2 = 1 / (N * 2 * Eb/N0)
// CP, pilot and preamble overhead are not counted (standard convention when
// comparing against theoretical BER curves).
double noiseVarianceForEbN0(const OfdmParams& p, double ebn0Linear);

// Theoretical QPSK BER (Gray mapping).
double qpskBerAwgn(double ebn0Linear);      // 0.5 * erfc(sqrt(Eb/N0))
double qpskBerRayleigh(double ebn0Linear);  // 0.5 * (1 - sqrt(g / (1 + g)))

inline double dbToLinear(double db) { return std::pow(10.0, db / 10.0); }

}  // namespace ofdm
