#include "ofdm/fft.hpp"

#include <cmath>
#include <map>
#include <stdexcept>
#include <utility>

namespace ofdm {

bool isPowerOfTwo(std::size_t n) { return n != 0 && (n & (n - 1)) == 0; }

namespace {

// Forward twiddles e^{-j2pi j/n}, j < n/2, computed once per FFT size.
// Each value is computed directly with cos/sin (no recurrence) for accuracy.
// Note: the cache is not thread-safe; add a mutex if FFTs are run from several threads.
const CVec& twiddles(std::size_t n) {
    static std::map<std::size_t, CVec> cache;
    auto it = cache.find(n);
    if (it == cache.end()) {
        CVec w(n / 2);
        for (std::size_t j = 0; j < n / 2; ++j) {
            const double angle = -2.0 * kPi * double(j) / double(n);
            w[j] = Complex(std::cos(angle), std::sin(angle));
        }
        it = cache.emplace(n, std::move(w)).first;
    }
    return it->second;
}

void fftInPlace(CVec& a, bool inverse) {
    const std::size_t n = a.size();
    if (!isPowerOfTwo(n)) throw std::invalid_argument("FFT: size must be a power of two");

    // 1) Bit-reversal permutation, so the butterflies can work in place.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }

    // 2) log2(N) stages of butterflies. Stage with block length `len` combines
    //    two DFTs of size len/2 into one DFT of size len.
    const CVec& tw = twiddles(n);
    for (std::size_t len = 2; len <= n; len <<= 1) {
        const std::size_t half = len / 2;
        const std::size_t step = n / len;  // e^{-j2pi j/len} = tw[j * n/len]
        for (std::size_t j = 0; j < half; ++j) {
            const Complex w = inverse ? std::conj(tw[j * step]) : tw[j * step];
            for (std::size_t i = 0; i < n; i += len) {
                const Complex u = a[i + j];
                const Complex v = a[i + j + half] * w;
                a[i + j] = u + v;
                a[i + j + half] = u - v;
            }
        }
    }

    if (inverse) {
        const double scale = 1.0 / double(n);
        for (Complex& v : a) v *= scale;
    }
}

CVec directTransform(const CVec& x, bool inverse) {
    const std::size_t n = x.size();
    const double sign = inverse ? 1.0 : -1.0;
    CVec out(n);
    for (std::size_t k = 0; k < n; ++k) {
        Complex sum(0.0, 0.0);
        for (std::size_t m = 0; m < n; ++m) {
            // (k*m) % n keeps the angle small for better accuracy at large N.
            const double angle = sign * 2.0 * kPi * double((k * m) % n) / double(n);
            sum += x[m] * Complex(std::cos(angle), std::sin(angle));
        }
        out[k] = inverse ? sum / double(n) : sum;
    }
    return out;
}

}  // namespace

CVec fft(CVec x) {
    fftInPlace(x, false);
    return x;
}

CVec ifft(CVec x) {
    fftInPlace(x, true);
    return x;
}

CVec dft(const CVec& x) { return directTransform(x, false); }
CVec idft(const CVec& x) { return directTransform(x, true); }

}  // namespace ofdm
