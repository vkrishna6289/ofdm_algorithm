#include "ofdm/pa.hpp"

#include <cmath>
#include <stdexcept>

#include "ofdm/metrics.hpp"

namespace ofdm {

CVec rappPa(const CVec& x, double saturation, double smoothness) {
    if (saturation <= 0.0) throw std::invalid_argument("rappPa: saturation must be > 0");
    if (smoothness <= 0.0) throw std::invalid_argument("rappPa: smoothness must be > 0");
    // Gain = (1 + u)^(-1/(2p)) with u = (|x|^2 / Asat^2)^p. It depends only on the
    // amplitude, so the phase is unchanged.
    // Fast path for p in {1, 2, 4, 8}: u^p by squaring and the (2p)-th root by
    // repeated sqrt, avoiding std::pow (very slow in MinGW's runtime).
    int sqrtCount = 0;  // 2p = 2^sqrtCount
    for (int s = 1; s <= 4; ++s)
        if (smoothness == double(1 << (s - 1))) sqrtCount = s;

    const double invSat2 = 1.0 / (saturation * saturation);
    CVec y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double r2 = mag2(x[i]) * invSat2;  // (|x| / Asat)^2
        double gain;
        if (sqrtCount > 0) {
            double u = r2;
            for (int s = 1; s < sqrtCount; ++s) u *= u;  // r2^p
            double root = 1.0 + u;
            for (int s = 0; s < sqrtCount; ++s) root = std::sqrt(root);
            gain = 1.0 / root;
        } else {
            gain = std::pow(1.0 + std::pow(r2, smoothness), -0.5 / smoothness);
        }
        y[i] = x[i] * gain;
    }
    return y;
}

double saturationForIbo(const CVec& x, double iboDb) {
    return std::sqrt(meanPower(x) * dbToLinear(iboDb));
}

}  // namespace ofdm
