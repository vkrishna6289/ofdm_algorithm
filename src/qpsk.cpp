#include "ofdm/qpsk.hpp"

#include <cmath>
#include <stdexcept>

namespace ofdm {

namespace {
const double kInvSqrt2 = 1.0 / std::sqrt(2.0);

void checkBit(int b) {
    if (b != 0 && b != 1) throw std::invalid_argument("QPSK: bit value must be 0 or 1");
}
}  // namespace

Complex qpskModulate(int b0, int b1) {
    checkBit(b0);
    checkBit(b1);
    const double i = (b0 == 0) ? 1.0 : -1.0;
    const double q = (b1 == 0) ? 1.0 : -1.0;
    return Complex(i * kInvSqrt2, q * kInvSqrt2);
}

CVec qpskModulate(const Bits& bits) {
    if (bits.size() % 2 != 0) throw std::invalid_argument("QPSK: bit count must be even");
    CVec symbols;
    symbols.reserve(bits.size() / 2);
    for (size_t i = 0; i < bits.size(); i += 2) symbols.push_back(qpskModulate(bits[i], bits[i + 1]));
    return symbols;
}

Bits qpskDemodulate(const CVec& symbols) {
    Bits bits;
    bits.reserve(symbols.size() * 2);
    for (const Complex& s : symbols) {
        bits.push_back(s.real() < 0.0 ? 1 : 0);
        bits.push_back(s.imag() < 0.0 ? 1 : 0);
    }
    return bits;
}

}  // namespace ofdm
