#pragma once

#include "ofdm/types.hpp"

namespace ofdm {

// Gray-mapped QPSK with unit average symbol energy:
//   s = (1/sqrt(2)) * (I + jQ),  bit 0 -> +1, bit 1 -> -1
// First bit of each pair drives I, second bit drives Q.
Complex qpskModulate(int b0, int b1);

// Modulates a bit vector (size must be even) into bits.size()/2 symbols.
CVec qpskModulate(const Bits& bits);

// Hard-decision demodulation: returns 2 bits per symbol.
Bits qpskDemodulate(const CVec& symbols);

}  // namespace ofdm
