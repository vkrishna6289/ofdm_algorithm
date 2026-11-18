#include "ofdm/channel.hpp"

#include <cmath>
#include <stdexcept>

#include "ofdm/fft.hpp"

namespace ofdm {

namespace {

// Complex value whose real and imaginary parts are independent N(0, 1),
// via Marsaglia's polar method on raw 53-bit uniforms.
// (std::normal_distribution is very slow with MinGW's libstdc++, because
// generate_canonical evaluates a long-double log on every call.)
Complex gaussianPair(std::mt19937_64& rng) {
    const double scale = 1.0 / 9007199254740992.0;  // 2^-53
    double u, v, s;
    do {
        u = 2.0 * double(rng() >> 11) * scale - 1.0;
        v = 2.0 * double(rng() >> 11) * scale - 1.0;
        s = u * u + v * v;
    } while (s >= 1.0 || s == 0.0);
    const double f = std::sqrt(-2.0 * std::log(s) / s);
    return Complex(u * f, v * f);
}

}  // namespace

CVec addAwgn(const CVec& x, double noiseVariance, std::mt19937_64& rng) {
    if (noiseVariance < 0.0) throw std::invalid_argument("addAwgn: negative noise variance");
    const double sigma = std::sqrt(noiseVariance / 2.0);  // per real dimension
    CVec y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = x[i] + sigma * gaussianPair(rng);
    return y;
}

CVec applyMultipath(const CVec& x, const CVec& taps) {
    CVec y(x.size(), Complex(0.0, 0.0));
    for (std::size_t n = 0; n < x.size(); ++n) {
        for (std::size_t l = 0; l < taps.size() && l <= n; ++l) y[n] += taps[l] * x[n - l];
    }
    return y;
}

CVec rayleighTaps(std::size_t numTaps, double decay, std::mt19937_64& rng) {
    if (numTaps == 0) throw std::invalid_argument("rayleighTaps: need at least one tap");
    if (decay <= 0.0) throw std::invalid_argument("rayleighTaps: decay must be positive");

    std::vector<double> power(numTaps);
    double total = 0.0;
    for (std::size_t l = 0; l < numTaps; ++l) {
        power[l] = std::exp(-double(l) / decay);
        total += power[l];
    }

    CVec taps(numTaps);
    for (std::size_t l = 0; l < numTaps; ++l) taps[l] = std::sqrt(power[l] / total / 2.0) * gaussianPair(rng);
    return taps;
}

CVec frequencyResponse(const CVec& taps, std::size_t nFft) {
    if (taps.size() > nFft) throw std::invalid_argument("frequencyResponse: more taps than nFft");
    CVec padded(nFft, Complex(0.0, 0.0));
    for (std::size_t l = 0; l < taps.size(); ++l) padded[l] = taps[l];
    return fft(padded);
}

}  // namespace ofdm
