// Unit tests for the OFDM library. Minimal self-contained framework so the
// project has no external dependencies. Run: ./ofdm_tests  (exit code 0 = pass)

#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "ofdm/channel.hpp"
#include "ofdm/fft.hpp"
#include "ofdm/metrics.hpp"
#include "ofdm/ofdm.hpp"
#include "ofdm/qpsk.hpp"

using namespace ofdm;

namespace {

int g_failures = 0;
int g_checks = 0;

#define CHECK(cond)                                                                       \
    do {                                                                                  \
        ++g_checks;                                                                       \
        if (!(cond)) {                                                                    \
            ++g_failures;                                                                 \
            std::cerr << "  FAIL " << __FILE__ << ":" << __LINE__ << "  " #cond << "\n"; \
        }                                                                                 \
    } while (0)

#define CHECK_THROWS(expr)                     \
    do {                                       \
        bool threw = false;                    \
        try {                                  \
            expr;                              \
        } catch (const std::exception&) {      \
            threw = true;                      \
        }                                      \
        CHECK(threw && "expected exception: " #expr); \
    } while (0)

double maxAbsDiff(const CVec& a, const CVec& b) {
    if (a.size() != b.size()) return 1e300;
    double m = 0.0;
    for (size_t i = 0; i < a.size(); ++i) m = std::max(m, std::abs(a[i] - b[i]));
    return m;
}

CVec randomVec(size_t n, std::mt19937_64& rng) {
    std::normal_distribution<double> g(0.0, 1.0);
    CVec v(n);
    for (auto& x : v) x = Complex(g(rng), g(rng));
    return v;
}

Bits randomBits(size_t n, std::mt19937_64& rng) {
    std::uniform_int_distribution<int> d(0, 1);
    Bits b(n);
    for (auto& x : b) x = d(rng);
    return b;
}

// ---- FFT ----------------------------------------------------------------

void testFftMatchesDft() {
    std::mt19937_64 rng(1);
    for (size_t n : {1u, 2u, 4u, 8u, 64u, 256u, 1024u}) {
        const CVec x = randomVec(n, rng);
        CHECK(maxAbsDiff(fft(x), dft(x)) < 1e-9 * double(n));
        CHECK(maxAbsDiff(ifft(x), idft(x)) < 1e-9);
    }
}

void testFftRoundTrip() {
    std::mt19937_64 rng(2);
    const CVec x = randomVec(512, rng);
    CHECK(maxAbsDiff(ifft(fft(x)), x) < 1e-12);
    CHECK(maxAbsDiff(fft(ifft(x)), x) < 1e-12);
}

void testFftKnownValues() {
    // Impulse -> flat spectrum
    CVec imp(8, 0.0);
    imp[0] = 1.0;
    CHECK(maxAbsDiff(fft(imp), CVec(8, 1.0)) < 1e-12);

    // Single complex exponential at bin 3 -> energy only in bin 3
    CVec tone(16);
    for (size_t n = 0; n < 16; ++n) tone[n] = std::polar(1.0, 2.0 * kPi * 3.0 * double(n) / 16.0);
    const CVec X = fft(tone);
    for (size_t k = 0; k < 16; ++k) CHECK(std::abs(X[k] - (k == 3 ? Complex(16.0) : Complex(0.0))) < 1e-9);

    // The original fft.cpp example input, checked against the direct DFT
    const CVec orig = {{2, 0}, {3, 4}, {4, 6}, {5, 6}};
    CHECK(maxAbsDiff(fft(orig), dft(orig)) < 1e-12);
}

void testParseval() {
    std::mt19937_64 rng(3);
    const CVec x = randomVec(128, rng);
    const CVec X = fft(x);
    double et = 0, ef = 0;
    for (auto& v : x) et += std::norm(v);
    for (auto& v : X) ef += std::norm(v);
    CHECK(std::abs(et - ef / 128.0) < 1e-9 * et);
}

void testFftRejectsBadSize() {
    CHECK_THROWS(fft(CVec(6)));
    CHECK_THROWS(fft(CVec()));
}

// ---- QPSK ---------------------------------------------------------------

void testQpsk() {
    for (int b0 = 0; b0 <= 1; ++b0)
        for (int b1 = 0; b1 <= 1; ++b1) {
            const Complex s = qpskModulate(b0, b1);
            CHECK(std::abs(std::norm(s) - 1.0) < 1e-12);  // unit energy
            const Bits r = qpskDemodulate({s});
            CHECK(r[0] == b0 && r[1] == b1);
        }
    // Gray mapping: adjacent constellation points differ in exactly one bit
    CHECK(qpskModulate(0, 0) == std::conj(qpskModulate(0, 1)));
    CHECK(qpskModulate(0, 0) == -std::conj(qpskModulate(1, 0)));

    CHECK_THROWS(qpskModulate(Bits{0, 1, 1}));
    CHECK_THROWS(qpskModulate(2, 0));

    std::mt19937_64 rng(4);
    const Bits bits = randomBits(1000, rng);
    CHECK(qpskDemodulate(qpskModulate(bits)) == bits);
}

// ---- OFDM building blocks ---------------------------------------------

void testParams() {
    const OfdmParams p = OfdmParams::wifiLike();
    p.validate();
    CHECK(p.numData() == 48);
    CHECK(p.pilotCarriers.size() == 4);
    CHECK(p.bin(-1) == 63 && p.bin(0) == 0 && p.bin(26) == 26 && p.bin(-26) == 38);

    OfdmParams bad = p;
    bad.dataCarriers.push_back(7);  // collides with a pilot
    CHECK_THROWS(bad.validate());
    bad = p;
    bad.nFft = 48;
    CHECK_THROWS(bad.validate());
}

void testCyclicPrefix() {
    std::mt19937_64 rng(5);
    const CVec s = randomVec(64, rng);
    const CVec withCp = addCyclicPrefix(s, 16);
    CHECK(withCp.size() == 80);
    for (size_t i = 0; i < 16; ++i) CHECK(withCp[i] == s[48 + i]);
    CHECK(removeCyclicPrefix(withCp, 16) == s);
}

void testSubcarrierMapping() {
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(6);
    const CVec data = qpskModulate(randomBits(p.bitsPerSymbol(), rng));
    const CVec bins = mapSubcarriers(p, data);
    CHECK(bins.size() == 64);
    CHECK(bins[0] == Complex(0.0));                           // DC null
    for (int k = 27; k <= 37; ++k) CHECK(bins[k] == Complex(0.0));  // guard band
    CHECK(bins[p.bin(21)] == Complex(-1.0));                  // pilot
    CHECK(demapSubcarriers(p, bins) == data);
    CHECK_THROWS(mapSubcarriers(p, CVec(47)));
}

// ---- End-to-end ---------------------------------------------------------

void testLoopbackIdealChannel() {
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(7);
    const Bits bits = randomBits(1001, rng);  // deliberately not a multiple of 96
    const CVec tx = transmitFrame(p, bits);
    CHECK(tx.size() == (1 + numDataSymbols(p, bits.size())) * 80);
    CHECK(receiveFrame(p, tx, bits.size()) == bits);
}

void testMultipathWithinCp() {
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(8);
    const Bits bits = randomBits(2000, rng);
    const CVec taps = rayleighTaps(8, 3.0, rng);  // 8 taps <= CP of 16
    const CVec rx = applyMultipath(transmitFrame(p, bits), taps);

    // LS estimate must match the true response on used carriers (no noise).
    const CVec hTrue = frequencyResponse(taps, p.nFft);
    const CVec hEst = estimateChannelLS(p, rx);
    for (int k : p.dataCarriers) CHECK(std::abs(hEst[p.bin(k)] - hTrue[p.bin(k)]) < 1e-9);

    CHECK(receiveFrame(p, rx, bits.size()) == bits);         // LS estimate
    CHECK(receiveFrame(p, rx, bits.size(), hTrue) == bits);  // perfect CSI
}

void testWithoutEqualiserFails() {
    // Sanity check that the equaliser is doing real work.
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(9);
    const Bits bits = randomBits(2000, rng);
    const CVec taps = rayleighTaps(8, 3.0, rng);
    const CVec rx = applyMultipath(transmitFrame(p, bits), taps);
    const Bits noEq = receiveFrame(p, rx, bits.size(), CVec(p.nFft, Complex(1.0)));
    CHECK(countBitErrors(noEq, bits) > 100);
}

// ---- Channel & metrics --------------------------------------------------

void testAwgnVariance() {
    std::mt19937_64 rng(10);
    const CVec n = addAwgn(CVec(200000, 0.0), 0.25, rng);
    double p = 0, re = 0;
    for (auto& v : n) {
        p += std::norm(v);
        re += v.real() * v.real();
    }
    p /= double(n.size());
    re /= double(n.size());
    CHECK(std::abs(p - 0.25) < 0.005);
    CHECK(std::abs(re - 0.125) < 0.003);  // half the power in each dimension
}

void testPapr() {
    CVec flat(64);
    for (size_t n = 0; n < 64; ++n) flat[n] = std::polar(1.0, 0.3 * double(n));
    CHECK(std::abs(paprDb(flat)) < 1e-9);  // constant envelope -> 0 dB

    CVec spike(64, 0.0);
    spike[5] = 1.0;
    CHECK(std::abs(paprDb(spike) - 10.0 * std::log10(64.0)) < 1e-9);

    // All subcarriers = 1 -> impulse in time -> PAPR = 10log10(N)
    CHECK(std::abs(paprDb(ifft(CVec(64, 1.0))) - 10.0 * std::log10(64.0)) < 1e-9);
}

void testOversampling() {
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(11);
    const CVec bins = mapSubcarriers(p, qpskModulate(randomBits(p.bitsPerSymbol(), rng)));
    const CVec x1 = ifft(bins);
    const CVec x4 = oversampledSymbol(bins, 4);
    CHECK(x4.size() == 256);
    // Every 4th oversampled sample equals the Nyquist-rate sample (up to the 1/L scale).
    for (size_t n = 0; n < 64; ++n) CHECK(std::abs(x4[4 * n] * 4.0 - x1[n]) < 1e-12);
    // Oversampling can only reveal higher peaks, never lower.
    CHECK(paprDb(x4) >= paprDb(x1) - 1e-9);
}

void testBerMatchesTheoryAwgn() {
    // Statistical check: simulated QPSK-OFDM BER at Eb/N0 = 4 dB should be
    // within 10% of 0.5*erfc(sqrt(Eb/N0)) ~= 1.25e-2.
    const OfdmParams p = OfdmParams::wifiLike();
    std::mt19937_64 rng(12);
    const double ebn0 = dbToLinear(4.0);
    const double nv = noiseVarianceForEbN0(p, ebn0);
    size_t errors = 0, total = 0;
    for (int f = 0; f < 400; ++f) {
        const Bits bits = randomBits(960, rng);
        const CVec rx = addAwgn(transmitFrame(p, bits), nv, rng);
        errors += countBitErrors(receiveFrame(p, rx, bits.size(), CVec(p.nFft, Complex(1.0))), bits);
        total += bits.size();
    }
    const double ber = double(errors) / double(total);
    const double theory = qpskBerAwgn(ebn0);
    std::cout << "    (BER @4dB sim=" << ber << " theory=" << theory << ")\n";
    CHECK(std::abs(ber - theory) / theory < 0.10);
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"fft matches dft", testFftMatchesDft},
        {"fft round trip", testFftRoundTrip},
        {"fft known values", testFftKnownValues},
        {"parseval", testParseval},
        {"fft rejects bad size", testFftRejectsBadSize},
        {"qpsk", testQpsk},
        {"ofdm params", testParams},
        {"cyclic prefix", testCyclicPrefix},
        {"subcarrier mapping", testSubcarrierMapping},
        {"loopback ideal channel", testLoopbackIdealChannel},
        {"multipath within CP", testMultipathWithinCp},
        {"no equaliser fails", testWithoutEqualiserFails},
        {"awgn variance", testAwgnVariance},
        {"papr", testPapr},
        {"oversampling", testOversampling},
        {"BER matches theory (AWGN)", testBerMatchesTheoryAwgn},
    };

    for (const auto& t : tests) {
        const int before = g_failures;
        try {
            t.second();
        } catch (const std::exception& e) {
            ++g_failures;
            std::cerr << "  EXCEPTION: " << e.what() << "\n";
        }
        std::cout << (g_failures == before ? "[PASS] " : "[FAIL] ") << t.first << "\n";
    }
    std::cout << "\n" << g_checks << " checks, " << g_failures << " failures\n";
    return g_failures == 0 ? 0 : 1;
}
