// PAPR CCDF of the baseline QPSK-OFDM transmitter — the reference curve that
// any PAPR-reduction scheme (SLM, PTS, clipping) must beat.
// Usage: papr_ccdf [numSymbols] [output.csv]   (defaults: 100000, papr_ccdf.csv)

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "ofdm/fft.hpp"
#include "ofdm/metrics.hpp"
#include "ofdm/ofdm.hpp"
#include "ofdm/qpsk.hpp"

using namespace ofdm;

namespace {

// Fraction of values strictly greater than threshold (values must be sorted).
double ccdfAt(const std::vector<double>& sorted, double threshold) {
    const auto it = std::upper_bound(sorted.begin(), sorted.end(), threshold);
    return double(sorted.end() - it) / double(sorted.size());
}

// PAPR value exceeded with probability `prob`.
double paprAtProbability(const std::vector<double>& sorted, double prob) {
    const std::size_t idx = static_cast<std::size_t>((1.0 - prob) * double(sorted.size()));
    return sorted[std::min(idx, sorted.size() - 1)];
}

}  // namespace

int main(int argc, char** argv) {
    const long numSymbols = argc > 1 ? std::atol(argv[1]) : 100000;
    const char* outPath = argc > 2 ? argv[2] : "papr_ccdf.csv";
    if (numSymbols <= 0) {
        std::cerr << "numSymbols must be positive\n";
        return 1;
    }

    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(7);
    std::uniform_int_distribution<int> coin(0, 1);

    std::vector<double> papr1, papr4;
    papr1.reserve(numSymbols);
    papr4.reserve(numSymbols);
    Bits bits(p.bitsPerSymbol());

    for (long s = 0; s < numSymbols; ++s) {
        for (auto& b : bits) b = coin(rng);
        const CVec bins = mapSubcarriers(p, qpskModulate(bits));
        papr1.push_back(paprDb(ifft(bins)));                 // Nyquist rate
        papr4.push_back(paprDb(oversampledSymbol(bins, 4)));  // 4x oversampled
    }
    std::sort(papr1.begin(), papr1.end());
    std::sort(papr4.begin(), papr4.end());

    // Classic approximation (independent Gaussian samples, Nyquist rate):
    // CCDF(z) = 1 - (1 - e^{-z})^N,  N = number of used subcarriers.
    const double nUsed = double(p.numData() + p.pilotCarriers.size());

    std::ofstream csv(outPath);
    if (!csv) {
        std::cerr << "cannot open " << outPath << "\n";
        return 1;
    }
    csv << "papr_db,ccdf_nyquist,ccdf_oversampled4,ccdf_theory_approx\n";
    for (double z = 0.0; z <= 13.0 + 1e-9; z += 0.1) {
        const double theory = 1.0 - std::pow(1.0 - std::exp(-dbToLinear(z)), nUsed);
        csv << z << ',' << ccdfAt(papr1, z) << ',' << ccdfAt(papr4, z) << ',' << theory << '\n';
    }

    std::cout << std::fixed << std::setprecision(2) << "Baseline QPSK-OFDM PAPR (" << numSymbols
              << " symbols, " << p.nFft << "-pt FFT, " << nUsed << " used subcarriers)\n\n"
              << "                     Nyquist   4x oversampled\n"
              << "mean PAPR [dB]     : " << std::setw(7) << [&] { double s = 0; for (double v : papr1) s += v; return s / papr1.size(); }()
              << "   " << std::setw(7) << [&] { double s = 0; for (double v : papr4) s += v; return s / papr4.size(); }() << "\n";
    for (double prob : {1e-1, 1e-2, 1e-3}) {
        if (prob * double(numSymbols) < 10) continue;  // not enough samples to estimate this tail
        std::cout << "PAPR @ CCDF=" << std::setw(5) << std::scientific << std::setprecision(0) << prob << std::fixed
                  << std::setprecision(2) << ": " << std::setw(7) << paprAtProbability(papr1, prob) << "   "
                  << std::setw(7) << paprAtProbability(papr4, prob) << "\n";
    }
    std::cout << "max PAPR [dB]      : " << std::setw(7) << papr1.back() << "   " << std::setw(7) << papr4.back()
              << "   (theoretical upper bound 10log10(" << nUsed << ") = " << 10.0 * std::log10(nUsed) << ")\n"
              << "\nWrote " << outPath << "\n";
    return 0;
}
