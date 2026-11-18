// Effect of a nonlinear power amplifier (Rapp model) on the baseline OFDM signal.
//
// Part 1 (noise-free): spectrum, ACLR and EVM vs input back-off (IBO).
// Part 2: BER vs Eb/N0 in AWGN for several IBOs.
//
// The signal is generated 4x oversampled so spectral regrowth falls into the
// extra guard band instead of aliasing on top of the signal.
//
// Usage: pa_sim [results_dir]   (default: current directory)
//   writes pa_psd.csv, pa_metrics.csv, pa_ber.csv

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include "ofdm/channel.hpp"
#include "ofdm/metrics.hpp"
#include "ofdm/ofdm.hpp"
#include "ofdm/pa.hpp"
#include "ofdm/qpsk.hpp"

using namespace ofdm;

namespace {

const std::size_t kOversampling = 4;
const std::size_t kWindow = 4;  // WOLA ramp length at the 1x rate (200 ns at 20 MHz)
const double kSmoothness = 2.0;  // Rapp p; typical SSPA
const double kLinear = std::numeric_limits<double>::infinity();  // "no PA" marker

// 802.11a-like channelisation: occupied band |k| <= 26, adjacent channel is
// the same width centred one channel spacing (64 subcarriers) away.
const int kMainLo = -26, kMainHi = 26;
const int kAdjLo = 64 - 26, kAdjHi = 64 + 26;

std::string label(double ibo) {
    if (std::isinf(ibo)) return "linear";
    std::ostringstream s;
    s << "ibo" << ibo;
    return s.str();
}

CVec amplify(const CVec& x, double ibo) {
    if (std::isinf(ibo)) return x;
    return rappPa(x, saturationForIbo(x, ibo), kSmoothness);
}

Bits randomBits(std::size_t n, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> coin(0, 1);
    Bits b(n);
    for (auto& v : b) v = coin(rng);
    return b;
}

// Least-squares complex gain alpha minimising |y - alpha x|^2 (Bussgang gain).
Complex linearGain(const CVec& x, const CVec& y) {
    Complex num(0.0, 0.0);
    double den = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        num += std::conj(x[i]) * y[i];
        den += mag2(x[i]);
    }
    return num / den;
}

}  // namespace

int main(int argc, char** argv) {
    const std::string dir = argc > 1 ? std::string(argv[1]) + "/" : "";
    OfdmParams base = OfdmParams::wifiLike();
    const OfdmParams rect = base.oversampled(kOversampling);  // no TX windowing, for reference
    base.windowLength = kWindow;
    const OfdmParams p = base.oversampled(kOversampling);
    std::mt19937_64 rng(99);

    // ---------------- Part 1: spectrum / ACLR / EVM ----------------
    // First entry: linear PA with rectangular (unwindowed) symbols, to show the
    // TX's own leakage floor. All other entries use the windowed transmitter.
    struct Config {
        std::string name;
        double ibo;
        bool windowed;
    };
    std::vector<Config> configs = {{"linear_rect", kLinear, false}, {"linear", kLinear, true}};
    for (double ibo : {12.0, 10.0, 8.0, 6.0, 4.0, 2.0, 0.0}) configs.push_back({label(ibo), ibo, true});

    const std::size_t nSym = 400;
    const Bits bits = randomBits(nSym * p.bitsPerSymbol(), rng);
    const CVec txWin = transmitFrame(p, bits);
    const CVec txRect = transmitFrame(rect, bits);
    const CVec ref = qpskModulate(bits);

    std::vector<std::vector<double>> psdDb;
    std::ofstream met(dir + "pa_metrics.csv");
    met << "config,ibo_db,obo_db,aclr_db,evm_pct,evm_db\n";

    std::cout << "Rapp PA, p = " << kSmoothness << ", " << kOversampling << "x oversampled, WOLA window "
              << kWindow << " samples, " << nSym << " OFDM symbols\n\n"
              << "   config    |  OBO   |  ACLR   |  EVM\n"
              << "             |  [dB]  |  [dB]   |  [%]    [dB]\n"
              << "-------------+--------+---------+----------------\n";
    for (const Config& c : configs) {
        const CVec& tx = c.windowed ? txWin : txRect;
        const double ibo = c.ibo;
        const CVec y = amplify(tx, ibo);

        const std::vector<double> psd = welchPsd(y, p.nFft);
        const double main = bandPower(psd, kMainLo, kMainHi);
        const double adj = std::max(bandPower(psd, kAdjLo, kAdjHi), bandPower(psd, -kAdjHi, -kAdjLo));
        const double aclr = 10.0 * std::log10(adj / main);

        const CVec eq = equalizeFrame(p, y, nSym, estimateChannelLS(p, y));
        const double evm = evmRms(eq, ref);

        // Output back-off: how far the mean output power sits below the PA's saturated output power.
        const double obo =
            std::isinf(ibo) ? kLinear : 10.0 * std::log10(std::pow(saturationForIbo(tx, ibo), 2) / meanPower(y));

        const double inbandMean = main / double(kMainHi - kMainLo + 1);
        std::vector<double> db(psd.size());
        for (std::size_t k = 0; k < psd.size(); ++k) db[k] = 10.0 * std::log10(std::max(psd[k] / inbandMean, 1e-30));
        psdDb.push_back(db);

        // EVM below -100 dB is numerical noise (no distortion); clamp for readability.
        const double evmDb = std::max(20.0 * std::log10(evm), -100.0);
        met << c.name << ',' << (std::isinf(ibo) ? -1 : ibo) << ',' << (std::isinf(obo) ? -1 : obo) << ',' << aclr
            << ',' << 100 * evm << ',' << evmDb << '\n';
        std::cout << std::fixed << std::setprecision(1) << std::setw(12) << c.name << " | ";
        if (std::isinf(obo))
            std::cout << "   -  ";
        else
            std::cout << std::setw(6) << obo;
        std::cout << " | " << std::setw(6) << aclr << "  | " << std::setw(6) << std::setprecision(2) << 100 * evm
                  << "  " << std::setw(6) << std::setprecision(1) << evmDb << "\n";
    }

    std::ofstream psdCsv(dir + "pa_psd.csv");
    psdCsv << "subcarrier";
    for (const Config& c : configs) psdCsv << ',' << c.name;
    psdCsv << '\n';
    const int half = static_cast<int>(p.nFft / 2);
    for (int k = -half; k < half; ++k) {
        psdCsv << k;
        for (const auto& db : psdDb) psdCsv << ',' << db[p.bin(k)];
        psdCsv << '\n';
    }

    // ---------------- Part 2: BER vs Eb/N0 ----------------
    // PA output is rescaled to the input's mean power so Eb/N0 refers to the
    // power actually transmitted, and the receiver is given the PA's linear
    // (Bussgang) gain exactly. Any BER loss is therefore due only to the
    // nonlinear distortion, not to power loss or channel-estimation error.
    const std::vector<double> berIbos = {kLinear, 8, 6, 4, 2, 0};
    const std::size_t frameSymbols = 20, minErrors = 300, minFrames = 50, maxBits = 2000000;

    std::ofstream ber(dir + "pa_ber.csv");
    ber << "ebn0_db,theory";
    for (double ibo : berIbos) ber << ',' << label(ibo);
    ber << '\n';

    std::cout << "\nBER vs Eb/N0 (AWGN)\n Eb/N0 |   theory  ";
    for (double ibo : berIbos) std::cout << std::setw(10) << label(ibo) << ' ';
    std::cout << "\n";

    // BER only falls as Eb/N0 rises, so once a configuration shows zero errors
    // in maxBits it is skipped (recorded as 0) at all higher Eb/N0 values.
    std::vector<bool> done(berIbos.size(), false);

    for (int db = 0; db <= 16; db += 2) {
        const double ebn0 = dbToLinear(db);
        const double nv = noiseVarianceForEbN0(p, ebn0);
        ber << db << ',' << qpskBerAwgn(ebn0);
        std::cout << std::setw(4) << db << " dB | " << std::scientific << std::setprecision(2) << qpskBerAwgn(ebn0)
                  << " ";

        for (std::size_t ci = 0; ci < berIbos.size(); ++ci) {
            const double ibo = berIbos[ci];
            std::size_t errors = 0, total = 0, frames = 0;
            while (!done[ci] && (errors < minErrors || frames < minFrames) && total < maxBits) {
                const Bits b = randomBits(frameSymbols * p.bitsPerSymbol(), rng);
                const CVec x = transmitFrame(p, b);
                CVec y = amplify(x, ibo);
                const double scale = std::sqrt(meanPower(x) / meanPower(y));
                for (auto& v : y) v *= scale;
                const CVec h(p.nFft, linearGain(x, y));
                errors += countBitErrors(receiveFrame(p, addAwgn(y, nv, rng), b.size(), h), b);
                total += b.size();
                ++frames;
            }
            if (errors == 0) done[ci] = true;
            const double v = total ? double(errors) / double(total) : 0.0;
            ber << ',' << v;
            std::cout << std::setw(10) << v << ' ';
        }
        ber << '\n';
        std::cout << std::fixed << "\n";
    }
    std::cout << "\n(0 = no errors in " << maxBits << " bits)\nWrote " << dir << "pa_psd.csv, " << dir
              << "pa_metrics.csv, " << dir << "pa_ber.csv\n";
    return 0;
}
