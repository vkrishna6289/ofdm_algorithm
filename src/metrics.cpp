#include "ofdm/metrics.hpp"

#include <cmath>
#include <stdexcept>

#include "ofdm/fft.hpp"

namespace ofdm {

double paprDb(const CVec& x) {
    if (x.empty()) throw std::invalid_argument("paprDb: empty signal");
    double peak = 0.0, sum = 0.0;
    for (const Complex& v : x) {
        const double pw = std::norm(v);  // |v|^2
        sum += pw;
        if (pw > peak) peak = pw;
    }
    if (sum == 0.0) throw std::invalid_argument("paprDb: zero-power signal");
    return 10.0 * std::log10(peak / (sum / double(x.size())));
}

CVec oversampledSymbol(const CVec& bins, std::size_t L) {
    const std::size_t n = bins.size();
    if (L == 0) throw std::invalid_argument("oversampledSymbol: L must be >= 1");
    if (!isPowerOfTwo(n) || !isPowerOfTwo(L)) throw std::invalid_argument("oversampledSymbol: sizes must be powers of two");

    // Positive frequencies stay at the start, negative frequencies move to the end,
    // zeros go in the middle (the new high-frequency region).
    CVec padded(n * L, Complex(0.0, 0.0));
    for (std::size_t k = 0; k < n / 2; ++k) padded[k] = bins[k];
    for (std::size_t k = n / 2; k < n; ++k) padded[n * L - n + k] = bins[k];
    return ifft(padded);
}

std::size_t countBitErrors(const Bits& a, const Bits& b) {
    if (a.size() != b.size()) throw std::invalid_argument("countBitErrors: size mismatch");
    std::size_t errors = 0;
    for (std::size_t i = 0; i < a.size(); ++i) errors += (a[i] != b[i]) ? 1 : 0;
    return errors;
}

double noiseVarianceForEbN0(const OfdmParams& p, double ebn0Linear) {
    if (ebn0Linear <= 0.0) throw std::invalid_argument("noiseVarianceForEbN0: Eb/N0 must be > 0");
    return 1.0 / (double(p.nFft) * 2.0 * ebn0Linear);
}

double qpskBerAwgn(double ebn0Linear) { return 0.5 * std::erfc(std::sqrt(ebn0Linear)); }

double qpskBerRayleigh(double ebn0Linear) {
    return 0.5 * (1.0 - std::sqrt(ebn0Linear / (1.0 + ebn0Linear)));
}

}  // namespace ofdm
