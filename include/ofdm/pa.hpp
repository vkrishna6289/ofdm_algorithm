#pragma once

#include "ofdm/types.hpp"

namespace ofdm {

// Rapp model of a solid-state power amplifier (AM/AM only, no AM/PM):
//   |y| = |x| / (1 + (|x| / Asat)^(2p))^(1/(2p)),   arg(y) = arg(x)
// Small signals pass with unit gain; the output amplitude saturates at Asat.
// p sets the sharpness of the knee: p ~ 2-3 for typical SSPAs, p -> inf is an
// ideal amplitude clipper.
CVec rappPa(const CVec& x, double saturation, double smoothness);

// Saturation amplitude that places the mean input power iboDb below the PA's
// input saturation power (input back-off, IBO):  Asat^2 = mean|x|^2 * 10^(IBO/10).
double saturationForIbo(const CVec& x, double iboDb);

}  // namespace ofdm
