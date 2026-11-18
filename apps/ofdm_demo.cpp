// Walks one frame through the full chain and prints what happens at each stage.

#include <iomanip>
#include <iostream>
#include <random>

#include "ofdm/channel.hpp"
#include "ofdm/fft.hpp"
#include "ofdm/metrics.hpp"
#include "ofdm/ofdm.hpp"
#include "ofdm/qpsk.hpp"

using namespace ofdm;

int main() {
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(2026);
    std::uniform_int_distribution<int> coin(0, 1);

    std::cout << "=== OFDM parameters ===\n"
              << "FFT size            : " << p.nFft << "\n"
              << "Cyclic prefix       : " << p.cpLength << " samples\n"
              << "Data subcarriers    : " << p.numData() << "\n"
              << "Pilot subcarriers   : " << p.pilotCarriers.size() << "\n"
              << "Bits per OFDM symbol: " << p.bitsPerSymbol() << " (QPSK)\n\n";

    // 1. Bits
    Bits bits(p.bitsPerSymbol() * 10);
    for (auto& b : bits) b = coin(rng);
    std::cout << "=== 1. Bits (first 16 of " << bits.size() << ") ===\n";
    for (size_t i = 0; i < 16; ++i) std::cout << bits[i];
    std::cout << "\n\n";

    // 2. QPSK
    const CVec sym = qpskModulate(Bits(bits.begin(), bits.begin() + 8));
    std::cout << std::fixed << std::setprecision(4) << "=== 2. QPSK (first 4 symbols) ===\n";
    for (const auto& s : sym) std::cout << "  " << s << "\n";
    std::cout << "\n";

    // 3-5. Mapping, IFFT, CP (done inside transmitFrame)
    std::vector<CVec> dataBins;
    const CVec tx = transmitFrame(p, bits, dataBins);
    std::cout << "=== 3-5. Mapping -> IFFT -> CP ===\n"
              << "Frame: 1 preamble + " << dataBins.size() << " data symbols = " << tx.size()
              << " complex samples\n";

    std::cout << "PAPR per data symbol (Nyquist / 4x oversampled), dB:\n";
    for (size_t s = 0; s < dataBins.size(); ++s)
        std::cout << "  symbol " << s << ": " << std::setw(6) << paprDb(ifft(dataBins[s])) << " / "
                  << std::setw(6) << paprDb(oversampledSymbol(dataBins[s], 4)) << "\n";
    std::cout << "\n";

    // 6. Channels + receiver
    std::cout << "=== 6. Channel + receiver ===\n";
    const double ebn0Db = 10.0;
    const double nv = noiseVarianceForEbN0(p, dbToLinear(ebn0Db));

    const Bits rxIdeal = receiveFrame(p, tx, bits.size());
    std::cout << "No channel           : " << countBitErrors(rxIdeal, bits) << " bit errors\n";

    const Bits rxAwgn = receiveFrame(p, addAwgn(tx, nv, rng), bits.size());
    std::cout << "AWGN @ " << ebn0Db << " dB Eb/N0 : " << countBitErrors(rxAwgn, bits) << " bit errors (theory BER "
              << std::scientific << qpskBerAwgn(dbToLinear(ebn0Db)) << std::fixed << ")\n";

    const CVec taps = rayleighTaps(6, 2.0, rng);
    const CVec rxMp = addAwgn(applyMultipath(tx, taps), nv, rng);
    std::cout << "Rayleigh 6-tap + AWGN: " << countBitErrors(receiveFrame(p, rxMp, bits.size()), bits)
              << " bit errors with LS equaliser, "
              << countBitErrors(receiveFrame(p, rxMp, bits.size(), CVec(p.nFft, Complex(1.0))), bits)
              << " without equaliser\n";
    return 0;
}
