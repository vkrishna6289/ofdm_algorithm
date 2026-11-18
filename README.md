# ofdm_algo — OFDM transmitter/receiver in C++ from scratch

A QPSK-OFDM physical layer built without DSP or communications libraries. It includes channel models, BER and PAPR
measurement, and serves as the baseline for low-PAPR research (SLM, PTS, clipping).

## Layout

```
include/ofdm/   public headers
  types.hpp     Complex, CVec, Bits
  qpsk.hpp      Gray-mapped unit-energy QPSK mod/demod
  fft.hpp       radix-2 FFT/IFFT + O(N^2) reference DFT
  ofdm.hpp      params, subcarrier mapping, CP, frame TX/RX, LS channel estimation
  channel.hpp   AWGN, FIR multipath, Rayleigh taps, frequency response
  metrics.hpp   PAPR, oversampling, BER helpers, theoretical BER
src/            implementations
apps/           ofdm_demo, ber_sim, papr_ccdf
tests/          unit tests (no external framework)
scripts/        plot_results.py
```

## System

| Item              | Value                                                        |
| ----------------- | ------------------------------------------------------------ |
| FFT size          | 64                                                           |
| Cyclic prefix     | 16 samples                                                   |
| Subcarriers       | 48 data, 4 pilots (±7, ±21), DC null, 11 guard (802.11a-style) |
| Modulation        | QPSK, Gray, unit energy                                      |
| Frame             | 1 preamble symbol (known QPSK on all 52 used carriers) + data symbols |
| Equalisation      | Zero-forcing using LS estimate from preamble, or perfect CSI |
| Synchronisation   | Ideal timing (not implemented yet)                           |

Conventions: `ifft` has the 1/N scale and `fft` has none. Eb/N0 is per data subcarrier, and CP, pilot and preamble
overhead are excluded (see `noiseVarianceForEbN0`).

## Build and run

```sh
cmake -S . -B build -G "MinGW Makefiles"   # omit -G on Linux
cmake --build build
./build/ofdm_tests                          # unit tests
./build/ofdm_demo                           # one frame through the whole chain
mkdir results
./build/ber_sim   results/ber.csv
./build/papr_ccdf 100000 results/papr_ccdf.csv
python scripts/plot_results.py results
```

## Baseline results

BER (QPSK-OFDM):
- AWGN matches `0.5·erfc(√(Eb/N0))`.
- Rayleigh with perfect CSI matches `0.5(1−√(γ/(1+γ)))`.
- With LS estimation from one preamble, the receiver is about 3 dB worse, because the estimate is as noisy as the data.

PAPR (100k symbols):

|                  | Nyquist | 4× oversampled |
| ---------------- | ------- | -------------- |
| mean             | 6.7 dB  | 7.2 dB         |
| PAPR @ CCDF 1e-3 | 10.3 dB | **10.6 dB**    |

**The 10.6 dB value is the number to beat.** Always measure PAPR with L ≥ 4 oversampling. Nyquist-rate samples miss
the true peaks.

## Roadmap

- [x] QPSK, FFT/IFFT (validated against the DFT)
- [x] Subcarrier mapping, CP, OFDM TX/RX
- [x] AWGN and Rayleigh multipath channels, LS channel estimation
- [x] BER vs theory, PAPR CCDF
- [ ] Spectrum / out-of-band emission measurement (PSD via Welch)
- [ ] Timing and CFO synchronisation, pilot-based phase tracking
- [ ] 16/64-QAM
- [ ] PAPR reduction: clipping and filtering, SLM, PTS, with BER and complexity trade-offs
- [ ] PA nonlinearity model (e.g. Rapp) so that PAPR reduction becomes visible in BER and spectrum
