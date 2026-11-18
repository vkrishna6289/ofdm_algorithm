#include "ofdm/ofdm.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "ofdm/fft.hpp"
#include "ofdm/qpsk.hpp"

namespace ofdm {

// ---- OfdmParams ---------------------------------------------------------

OfdmParams OfdmParams::wifiLike() {
    OfdmParams p;
    p.nFft = 64;
    p.cpLength = 16;
    p.pilotCarriers = {-21, -7, 7, 21};
    p.pilotValues = {1.0, 1.0, 1.0, -1.0};
    for (int k = -26; k <= 26; ++k) {
        if (k == 0) continue;  // DC null
        if (std::find(p.pilotCarriers.begin(), p.pilotCarriers.end(), k) != p.pilotCarriers.end()) continue;
        p.dataCarriers.push_back(k);
    }
    return p;
}

std::size_t OfdmParams::bin(int k) const {
    const long n = static_cast<long>(nFft);
    return static_cast<std::size_t>(((k % n) + n) % n);
}

OfdmParams OfdmParams::oversampled(std::size_t L) const {
    if (!isPowerOfTwo(L)) throw std::invalid_argument("OfdmParams::oversampled: L must be a power of two");
    OfdmParams q = *this;
    q.nFft = nFft * L;
    q.cpLength = cpLength * L;
    q.windowLength = windowLength * L;
    return q;
}

void OfdmParams::validate() const {
    if (!isPowerOfTwo(nFft)) throw std::invalid_argument("OfdmParams: nFft must be a power of two");
    if (cpLength >= nFft) throw std::invalid_argument("OfdmParams: cpLength must be < nFft");
    if (windowLength > cpLength) throw std::invalid_argument("OfdmParams: windowLength must be <= cpLength");
    if (dataCarriers.empty()) throw std::invalid_argument("OfdmParams: no data carriers");
    if (pilotValues.size() != pilotCarriers.size())
        throw std::invalid_argument("OfdmParams: pilotValues/pilotCarriers size mismatch");

    const int half = static_cast<int>(nFft / 2);
    std::vector<bool> used(nFft, false);
    auto claim = [&](int k) {
        if (k < -half || k >= half) throw std::invalid_argument("OfdmParams: carrier index out of range");
        const std::size_t b = bin(k);
        if (used[b]) throw std::invalid_argument("OfdmParams: carrier used twice");
        used[b] = true;
    };
    for (int k : dataCarriers) claim(k);
    for (int k : pilotCarriers) claim(k);
}

// ---- Symbol-level -------------------------------------------------------

CVec mapSubcarriers(const OfdmParams& p, const CVec& dataSymbols) {
    if (dataSymbols.size() != p.numData())
        throw std::invalid_argument("mapSubcarriers: need exactly numData() symbols");
    CVec bins(p.nFft, Complex(0.0, 0.0));
    for (std::size_t i = 0; i < p.dataCarriers.size(); ++i) bins[p.bin(p.dataCarriers[i])] = dataSymbols[i];
    for (std::size_t i = 0; i < p.pilotCarriers.size(); ++i) bins[p.bin(p.pilotCarriers[i])] = p.pilotValues[i];
    return bins;
}

CVec demapSubcarriers(const OfdmParams& p, const CVec& bins) {
    if (bins.size() != p.nFft) throw std::invalid_argument("demapSubcarriers: need nFft bins");
    CVec data;
    data.reserve(p.numData());
    for (int k : p.dataCarriers) data.push_back(bins[p.bin(k)]);
    return data;
}

CVec addCyclicPrefix(const CVec& symbol, std::size_t cpLength) {
    if (cpLength > symbol.size()) throw std::invalid_argument("addCyclicPrefix: CP longer than symbol");
    CVec out;
    out.reserve(symbol.size() + cpLength);
    out.insert(out.end(), symbol.end() - static_cast<long>(cpLength), symbol.end());
    out.insert(out.end(), symbol.begin(), symbol.end());
    return out;
}

CVec removeCyclicPrefix(const CVec& symbolWithCp, std::size_t cpLength) {
    if (cpLength > symbolWithCp.size()) throw std::invalid_argument("removeCyclicPrefix: input too short");
    return CVec(symbolWithCp.begin() + static_cast<long>(cpLength), symbolWithCp.end());
}

CVec modulateSymbol(const OfdmParams& p, const CVec& bins) {
    if (bins.size() != p.nFft) throw std::invalid_argument("modulateSymbol: need nFft bins");
    return addCyclicPrefix(ifft(bins), p.cpLength);
}

CVec demodulateSymbol(const OfdmParams& p, const CVec& samples) {
    if (samples.size() != p.symbolLength()) throw std::invalid_argument("demodulateSymbol: wrong length");
    return fft(removeCyclicPrefix(samples, p.cpLength));
}

// ---- Frame-level --------------------------------------------------------

CVec preambleBins(const OfdmParams& p) {
    // Deterministic pseudo-random QPSK on all used carriers. mt19937's output
    // sequence is fixed by the C++ standard, so TX and RX always agree.
    std::mt19937 rng(0x0FD3u);
    CVec bins(p.nFft, Complex(0.0, 0.0));
    auto fill = [&](int k) {
        const int b0 = static_cast<int>((rng() >> 31) & 1u);
        const int b1 = static_cast<int>((rng() >> 31) & 1u);
        bins[p.bin(k)] = qpskModulate(b0, b1);
    };
    for (int k : p.dataCarriers) fill(k);
    for (int k : p.pilotCarriers) fill(k);
    return bins;
}

std::size_t numDataSymbols(const OfdmParams& p, std::size_t numBits) {
    const std::size_t bps = p.bitsPerSymbol();
    return (numBits + bps - 1) / bps;
}

namespace {

// Adds OFDM symbol number `index` (0 = preamble) into the frame buffer, with
// optional WOLA windowing. Symbol `index` starts at index * symbolLength().
void addSymbolToFrame(const OfdmParams& p, const CVec& bins, std::size_t index, CVec& frame) {
    const std::size_t w = p.windowLength;
    const CVec body = ifft(bins);
    CVec ext = addCyclicPrefix(body, p.cpLength);
    ext.insert(ext.end(), body.begin(), body.begin() + static_cast<long>(w));  // cyclic suffix

    // Raised-cosine ramps; rampUp[n] + rampDown[n] = 1, so the overlapping
    // tail of one symbol and head of the next sum to a smooth transition.
    for (std::size_t n = 0; n < w; ++n) {
        const double up = 0.5 * (1.0 - std::cos(kPi * (double(n) + 0.5) / double(w)));
        ext[n] *= up;
        ext[ext.size() - w + n] *= 1.0 - up;
    }

    const std::size_t start = index * p.symbolLength();
    for (std::size_t n = 0; n < ext.size(); ++n) frame[start + n] += ext[n];
}

}  // namespace

CVec transmitFrame(const OfdmParams& p, const Bits& bits, std::vector<CVec>& dataBinsOut) {
    p.validate();
    const std::size_t bps = p.bitsPerSymbol();
    const std::size_t nSym = numDataSymbols(p, bits.size());

    CVec frame((nSym + 1) * p.symbolLength() + p.windowLength, Complex(0.0, 0.0));
    addSymbolToFrame(p, preambleBins(p), 0, frame);

    dataBinsOut.clear();
    dataBinsOut.reserve(nSym);
    for (std::size_t s = 0; s < nSym; ++s) {
        Bits chunk(bps, 0);  // zero padding for the last symbol
        const std::size_t start = s * bps;
        const std::size_t count = std::min(bps, bits.size() - start);
        std::copy(bits.begin() + static_cast<long>(start), bits.begin() + static_cast<long>(start + count),
                  chunk.begin());

        dataBinsOut.push_back(mapSubcarriers(p, qpskModulate(chunk)));
        addSymbolToFrame(p, dataBinsOut.back(), s + 1, frame);
    }
    return frame;
}

CVec transmitFrame(const OfdmParams& p, const Bits& bits) {
    std::vector<CVec> unused;
    return transmitFrame(p, bits, unused);
}

namespace {
CVec symbolAt(const OfdmParams& p, const CVec& frame, std::size_t index) {
    const std::size_t len = p.symbolLength();
    const std::size_t start = index * len;
    if (start + len > frame.size()) throw std::invalid_argument("receiver: frame too short");
    return CVec(frame.begin() + static_cast<long>(start), frame.begin() + static_cast<long>(start + len));
}
}  // namespace

CVec estimateChannelLS(const OfdmParams& p, const CVec& frameSamples) {
    p.validate();
    const CVec y = demodulateSymbol(p, symbolAt(p, frameSamples, 0));
    const CVec x = preambleBins(p);
    CVec h(p.nFft, Complex(0.0, 0.0));
    for (int k : p.dataCarriers) h[p.bin(k)] = y[p.bin(k)] / x[p.bin(k)];
    for (int k : p.pilotCarriers) h[p.bin(k)] = y[p.bin(k)] / x[p.bin(k)];
    return h;
}

CVec equalizeFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numSymbols, const CVec& channel) {
    p.validate();
    if (channel.size() != p.nFft) throw std::invalid_argument("equalizeFrame: channel must have nFft bins");

    CVec out;
    out.reserve(numSymbols * p.numData());
    for (std::size_t s = 0; s < numSymbols; ++s) {
        CVec bins = demodulateSymbol(p, symbolAt(p, frameSamples, s + 1));
        for (int k : p.dataCarriers) bins[p.bin(k)] /= channel[p.bin(k)];  // zero-forcing equaliser
        const CVec data = demapSubcarriers(p, bins);
        out.insert(out.end(), data.begin(), data.end());
    }
    return out;
}

Bits receiveFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numBits, const CVec& channel) {
    Bits out = qpskDemodulate(equalizeFrame(p, frameSamples, numDataSymbols(p, numBits), channel));
    out.resize(numBits);
    return out;
}

Bits receiveFrame(const OfdmParams& p, const CVec& frameSamples, std::size_t numBits) {
    return receiveFrame(p, frameSamples, numBits, estimateChannelLS(p, frameSamples));
}

}  // namespace ofdm
