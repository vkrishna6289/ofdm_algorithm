#pragma once

#include <cstddef>
#include <vector>

#include "ofdm/types.hpp"

namespace ofdm {

// Describes the OFDM numerology and how subcarriers are used.
// Subcarriers are addressed by their logical index k in [-N/2, N/2),
// where k = 0 is DC. bin(k) converts to the FFT bin index in [0, N).
struct OfdmParams {
    std::size_t nFft = 64;
    std::size_t cpLength = 16;
    // Transmit windowing (WOLA). Each symbol is extended by a cyclic suffix of
    // windowLength samples, its first/last windowLength samples get raised-cosine
    // ramps, and consecutive symbols overlap-add over that region. This smooths
    // the symbol boundaries and strongly reduces out-of-band emission.
    // The overlap lies inside the CP, so the receiver is unchanged, but the
    // multipath delay spread the CP can absorb shrinks to cpLength - windowLength.
    // 0 = no windowing (rectangular symbols). Must be <= cpLength.
    std::size_t windowLength = 0;
    std::vector<int> dataCarriers;
    std::vector<int> pilotCarriers;
    CVec pilotValues;  // one known value per pilot carrier

    // 802.11a-style layout: N=64, CP=16, 48 data, 4 pilots at +-7, +-21,
    // DC null, guard bands on k = -32..-27 and 27..31 (52 used subcarriers).
    static OfdmParams wifiLike();

    std::size_t numData() const { return dataCarriers.size(); }
    std::size_t bitsPerSymbol() const { return 2 * numData(); }  // QPSK
    std::size_t symbolLength() const { return nFft + cpLength; }  // samples incl. CP
    std::size_t bin(int k) const;

    // Same system sampled L times faster: nFft and cpLength scale by L, the
    // occupied subcarriers (and the symbol duration) stay the same, and the
    // extra bins become wide guard bands. Use this whenever the signal passes a
    // nonlinearity (e.g. a PA) or its spectrum is measured; at the Nyquist rate
    // out-of-band regrowth would alias back on top of the signal.
    OfdmParams oversampled(std::size_t L) const;

    // Throws std::invalid_argument if the configuration is inconsistent.
    void validate() const;
};

// ---- Single OFDM symbol building blocks --------------------------------

// Places data symbols and pilots into an N-bin frequency-domain vector.
// Unused bins (DC, guards) are zero.
CVec mapSubcarriers(const OfdmParams& p, const CVec& dataSymbols);

// Extracts the data-carrier values from an N-bin frequency-domain vector.
CVec demapSubcarriers(const OfdmParams& p, const CVec& bins);

CVec addCyclicPrefix(const CVec& symbol, std::size_t cpLength);
CVec removeCyclicPrefix(const CVec& symbolWithCp, std::size_t cpLength);

// bins -> IFFT -> add CP. Output has p.symbolLength() samples.
CVec modulateSymbol(const OfdmParams& p, const CVec& bins);

// remove CP -> FFT -> bins. Input must have p.symbolLength() samples.
CVec demodulateSymbol(const OfdmParams& p, const CVec& samples);

// ---- Frame-level transmitter / receiver ---------------------------------
//
// Frame = [preamble symbol][data symbol 0]...[data symbol M-1]
// The preamble carries a known QPSK sequence on every used subcarrier and is
// used by the receiver for least-squares channel estimation.
// Timing synchronisation is assumed ideal (the receiver knows where the frame starts).

// Known preamble in the frequency domain (N bins).
CVec preambleBins(const OfdmParams& p);

// Number of data OFDM symbols needed for numBits (last symbol zero-padded).
std::size_t numDataSymbols(const OfdmParams& p, std::size_t numBits);

// bits -> QPSK -> subcarrier mapping -> IFFT -> CP (-> windowing), with preamble in front.
// Output length: (1 + numDataSymbols) * symbolLength() + windowLength.
CVec transmitFrame(const OfdmParams& p, const Bits& bits);

// Same, but also returns the frequency-domain bins of every data symbol
// (useful for PAPR work, which operates on individual OFDM symbols).
CVec transmitFrame(const OfdmParams& p, const Bits& bits, std::vector<CVec>& dataBinsOut);

// Least-squares channel estimate H[k] = Y[k] / P[k] from the preamble.
// Returns N bins; unused bins are zero.
CVec estimateChannelLS(const OfdmParams& p, const CVec& frameSamples);

// FFT + zero-forcing equalisation of the first numSymbols data symbols.
// Returns the equalised data-carrier values, numSymbols * numData() of them
// (before the QPSK decision; useful for EVM and constellation plots).
CVec equalizeFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numSymbols, const CVec& channel);

// Demodulates a frame using the supplied channel frequency response H (N bins)
// for zero-forcing equalisation, then returns the first numBits bits.
Bits receiveFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numBits, const CVec& channel);

// Same, estimating the channel from the preamble (LS).
Bits receiveFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numBits);

}  // namespace ofdm
