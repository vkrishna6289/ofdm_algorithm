#include "ofdm/channel.hpp"

#include <cmath>
#include <stdexcept>

#include "ofdm/fft.hpp"

namespace ofdm {

CVec addAwgn(const CVec& x, double noiseVariance, std::mt19937_64& rng) {
    if (noiseVariance < 0.0) throw std::invalid_argument("addAwgn: negative noise variance");
    std::normal_distribution<double> gauss(0.0, std::sqrt(noiseVariance / 2.0));
    CVec y(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) y[i] = x[i] + Complex(gauss(rng), gauss(rng));
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

    std::normal_distribution<double> gauss(0.0, 1.0);
    CVec taps(numTaps);
    for (std::size_t l = 0; l < numTaps; ++l) {
        const double sigma = std::sqrt(power[l] / total / 2.0);
        taps[l] = Complex(sigma * gauss(rng), sigma * gauss(rng));
    }
    return taps;
}

CVec frequencyResponse(const CVec& taps, std::size_t nFft) {
    if (taps.size() > nFft) throw std::invalid_argument("frequencyResponse: more taps than nFft");
    CVec padded(nFft, Complex(0.0, 0.0));
    for (std::size_t l = 0; l < taps.size(); ++l) padded[l] = taps[l];
    return fft(padded);
}

}  // namespace ofdm
