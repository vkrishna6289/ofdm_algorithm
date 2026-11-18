# ofdm_algo — OFDM transmitter/receiver in C++ from scratch

A QPSK-OFDM physical layer built without DSP or communications libraries. It includes channel models, a nonlinear PA
model, and BER, PAPR, spectrum and EVM measurement. It serves as the baseline for low-PAPR research (SLM, PTS,
clipping).

## Layout

```
include/ofdm/   public headers
  types.hpp     Complex, CVec, Bits, mag2
  qpsk.hpp      Gray-mapped unit-energy QPSK mod/demod
  fft.hpp       radix-2 FFT/IFFT (cached twiddles) + O(N^2) reference DFT
  ofdm.hpp      params, oversampling, subcarrier mapping, CP, WOLA windowing,
                frame TX/RX, LS channel estimation, ZF equalisation
  channel.hpp   AWGN, FIR multipath, Rayleigh taps, frequency response
  pa.hpp        Rapp power-amplifier model, input back-off
  metrics.hpp   PAPR, Welch PSD, band power / ACLR, EVM, BER helpers, theoretical BER
src/            implementations
apps/           ofdm_demo, ber_sim, papr_ccdf, pa_sim
tests/          unit tests (no external framework)
scripts/        plot_results.py
```

## System

| Item              | Value                                                        |
| ----------------- | ------------------------------------------------------------ |
| FFT size          | 64 (×L when oversampled with `OfdmParams::oversampled(L)`)   |
| Cyclic prefix     | 16 samples                                                   |
| Subcarriers       | 48 data, 4 pilots (±7, ±21), DC null, 11 guard (802.11a-style) |
| Modulation        | QPSK, Gray, unit energy                                      |
| TX windowing      | Optional WOLA raised-cosine (`windowLength`, ≤ CP; `pa_sim` uses 4) |
| Frame             | 1 preamble symbol (known QPSK on all 52 used carriers) + data symbols |
| Equalisation      | Zero-forcing using LS estimate from preamble, or perfect CSI |
| Synchronisation   | Ideal timing (not implemented yet)                           |

Conventions:
- `ifft` has the 1/N scale and `fft` has none.
- Eb/N0 is per data subcarrier, and CP, pilot and preamble overhead are excluded (see `noiseVarianceForEbN0`).
- Anything that passes through a nonlinearity or has its spectrum measured must use the oversampled system (L = 4);
  otherwise spectral regrowth aliases back onto the signal.

## Build and run

```sh
cmake -S . -B build -G "MinGW Makefiles"   # omit -G on Linux
cmake --build build
./build/ofdm_tests                          # unit tests (<1 s)
./build/ofdm_demo                           # one frame through the whole chain
mkdir results
./build/ber_sim   results/ber.csv           # ~15 s
./build/papr_ccdf 100000 results/papr_ccdf.csv  # ~1 s
./build/pa_sim    results                   # ~25 s
python scripts/plot_results.py results
```

Performance notes for MinGW/libstdc++:
- `-fcx-limited-range` stops complex multiplication from going through the slow `__muldc3` routine.
- `mag2()` replaces `std::norm`, which libstdc++ computes through `hypot`.
- The Gaussian generator replaces `std::normal_distribution`, which recomputes an extended-precision `log` on every
  call.

Together these make the simulations 10–15× faster.

## Baseline results

### BER (QPSK-OFDM)
- AWGN matches `0.5·erfc(√(Eb/N0))`.
- Rayleigh with perfect CSI matches `0.5(1−√(γ/(1+γ)))`.
- With LS estimation from one preamble, the receiver is about 3 dB worse, because the estimate is as noisy as the data.

### PAPR (100k symbols)

|                  | Nyquist | 4× oversampled |
| ---------------- | ------- | -------------- |
| mean             | 6.7 dB  | 7.2 dB         |
| PAPR @ CCDF 1e-3 | 10.3 dB | **10.6 dB**    |

**The 10.6 dB value is the number to beat.** Always measure PAPR with L ≥ 4 oversampling. Nyquist-rate samples miss
the true peaks.

### Power amplifier (Rapp, p = 2, 4× oversampled, WOLA window 4)

| Config                | OBO [dB] | ACLR [dB] | EVM [%] |
| --------------------- | -------- | --------- | ------- |
| linear, no TX window  | –        | −29.9     | 0       |
| linear                | –        | −48.1     | 0       |
| IBO 12 dB             | 12.0     | −46.7     | 0.7     |
| IBO 10 dB             | 10.1     | −43.5     | 1.6     |
| IBO 8 dB              | 8.3      | −38.7     | 3.4     |
| IBO 6 dB              | 6.5      | −33.8     | 6.5     |
| IBO 4 dB              | 5.0      | −29.6     | 10.9    |
| IBO 2 dB              | 3.7      | −26.1     | 16.5    |
| IBO 0 dB              | 2.7      | −23.4     | 22.8    |

ACLR is the power in the adjacent 802.11 channel (64 subcarriers away) relative to the main channel.

- Without TX windowing, the abrupt symbol boundaries alone limit ACLR to about −30 dB, which hides most of the PA's
  regrowth. Windowing is a prerequisite for meaningful spectral comparisons.
- Below about 6 dB IBO, distortion rises quickly. For reference, 802.11a allows −10 dB EVM for QPSK rate 1/2.
- BER (with the PA output renormalised and its linear gain known) is unaffected at 8 dB IBO and loses about 1 dB at
  1e-4 for 4 dB IBO. At 0 dB IBO it develops a distortion-limited tail.

**Why this matters for PAPR reduction:** the PA has to be backed off by roughly the PAPR to stay linear, and back-off
wastes power (OBO ≈ 8 dB means the PA's average output is about 1/6 of its capability). A scheme that cuts PAPR by
X dB lets the PA run with roughly X dB less back-off for the same ACLR and EVM. That trade-off is how SLM, PTS and
clipping should be compared.

## Roadmap

- [x] QPSK, FFT/IFFT (validated against the DFT)
- [x] Subcarrier mapping, CP, OFDM TX/RX
- [x] AWGN and Rayleigh multipath channels, LS channel estimation
- [x] BER vs theory, PAPR CCDF
- [x] Oversampled system, WOLA TX windowing
- [x] Rapp PA model; spectrum (Welch PSD), ACLR, EVM, BER vs back-off
- [ ] PAPR reduction: clipping and filtering, SLM, PTS; compare at equal ACLR/EVM, plus BER and complexity
- [ ] Timing and CFO synchronisation, pilot-based phase tracking
- [ ] 16/64-QAM (much more sensitive to PA distortion than QPSK)
