// BER vs Eb/N0 for QPSK-OFDM in AWGN and Rayleigh multipath, compared with theory.
// Usage: ber_sim [output.csv]   (default: ber.csv)

#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

#include "ofdm/channel.hpp"
#include "ofdm/metrics.hpp"
#include "ofdm/ofdm.hpp"

using namespace ofdm;

namespace {

struct Counter {
    std::size_t errors = 0, total = 0;
    double ber() const { return total ? double(errors) / double(total) : 0.0; }
};

const std::size_t kMinErrors = 500;
// Block fading draws one channel per frame, so many frames are needed to
// average over channel realisations even when errors accumulate quickly.
const std::size_t kMinFrames = 500;
const std::size_t kMaxBits = 5000000;
const std::size_t kSymbolsPerFrame = 20;
const std::size_t kChannelTaps = 6;   // <= CP length (16), so no ISI
const double kTapDecay = 2.0;

}  // namespace

int main(int argc, char** argv) {
    const char* outPath = argc > 1 ? argv[1] : "ber.csv";
    const OfdmParams p = OfdmParams::wifiLike();
    const std::size_t frameBits = p.bitsPerSymbol() * kSymbolsPerFrame;
    const CVec flat(p.nFft, Complex(1.0));

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int> coin(0, 1);

    std::ofstream csv(outPath);
    if (!csv) {
        std::cerr << "cannot open " << outPath << "\n";
        return 1;
    }
    csv << "ebn0_db,awgn_sim,awgn_theory,rayleigh_perfect_csi,rayleigh_ls,rayleigh_theory\n";

    std::cout << " Eb/N0 |  AWGN sim   AWGN theory | Rayl. perfCSI  Rayl. LS   Rayl. theory\n"
              << "-------+-------------------------+---------------------------------------\n";

    for (int db = 0; db <= 20; db += 2) {
        const double ebn0 = dbToLinear(db);
        const double nv = noiseVarianceForEbN0(p, ebn0);
        Counter awgn, perfect, ls;
        std::size_t frames = 0;

        while ((frames < kMinFrames || awgn.errors < kMinErrors || perfect.errors < kMinErrors ||
                ls.errors < kMinErrors) &&
               awgn.total < kMaxBits) {
            ++frames;
            Bits bits(frameBits);
            for (auto& b : bits) b = coin(rng);
            const CVec tx = transmitFrame(p, bits);

            // AWGN only
            const Bits rA = receiveFrame(p, addAwgn(tx, nv, rng), bits.size(), flat);
            awgn.errors += countBitErrors(rA, bits);
            awgn.total += bits.size();

            // Rayleigh multipath (new channel per frame = block fading) + AWGN
            const CVec taps = rayleighTaps(kChannelTaps, kTapDecay, rng);
            const CVec rx = addAwgn(applyMultipath(tx, taps), nv, rng);
            perfect.errors += countBitErrors(receiveFrame(p, rx, bits.size(), frequencyResponse(taps, p.nFft)), bits);
            perfect.total += bits.size();
            ls.errors += countBitErrors(receiveFrame(p, rx, bits.size()), bits);
            ls.total += bits.size();
        }

        csv << db << ',' << awgn.ber() << ',' << qpskBerAwgn(ebn0) << ',' << perfect.ber() << ',' << ls.ber() << ','
            << qpskBerRayleigh(ebn0) << '\n';
        std::cout << std::setw(4) << db << " dB | " << std::scientific << std::setprecision(3) << std::setw(10)
                  << awgn.ber() << "  " << std::setw(10) << qpskBerAwgn(ebn0) << "  | " << std::setw(10)
                  << perfect.ber() << "  " << std::setw(10) << ls.ber() << "  " << std::setw(10)
                  << qpskBerRayleigh(ebn0) << std::fixed << "\n";
    }
    std::cout << "\n(AWGN sim shows 0 once BER drops below ~1/" << kMaxBits << ")\nWrote " << outPath << "\n";
    return 0;
}
