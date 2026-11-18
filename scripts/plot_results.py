"""Plot BER and PAPR-CCDF curves from the CSVs written by ber_sim and papr_ccdf.

Usage: python scripts/plot_results.py [results_dir]   (default: results)
"""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return {k: [float(r[k]) for r in rows] for k in rows[0]}


def masked(xs, ys):
    """Drop zero values, which cannot be shown on a log axis."""
    pairs = [(x, y) for x, y in zip(xs, ys) if y > 0]
    return [p[0] for p in pairs], [p[1] for p in pairs]


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "results")

    ber = load(out / "ber.csv")
    x = ber["ebn0_db"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(*masked(x, ber["awgn_theory"]), "k-", label="AWGN theory")
    ax.semilogy(*masked(x, ber["awgn_sim"]), "o", label="AWGN sim")
    ax.semilogy(*masked(x, ber["rayleigh_theory"]), "k--", label="Rayleigh theory")
    ax.semilogy(*masked(x, ber["rayleigh_perfect_csi"]), "s", label="Rayleigh sim, perfect CSI")
    ax.semilogy(*masked(x, ber["rayleigh_ls"]), "^", label="Rayleigh sim, LS estimate")
    ax.set(xlabel="Eb/N0 [dB]", ylabel="BER", title="QPSK-OFDM BER (N=64, CP=16)", ylim=(1e-6, 0.5))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "ber.png", dpi=150)

    papr = load(out / "papr_ccdf.csv")
    z = papr["papr_db"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(*masked(z, papr["ccdf_theory_approx"]), "k--", label="Approx. theory (Nyquist)")
    ax.semilogy(*masked(z, papr["ccdf_nyquist"]), "-", label="Sim, Nyquist rate (L=1)")
    ax.semilogy(*masked(z, papr["ccdf_oversampled4"]), "-", label="Sim, 4x oversampled (L=4)")
    ax.set(xlabel="PAPR threshold z [dB]", ylabel="Pr(PAPR > z)", title="Baseline OFDM PAPR CCDF",
           ylim=(1e-4, 1), xlim=(4, 13))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "papr_ccdf.png", dpi=150)
    print(f"Wrote {out / 'ber.png'} and {out / 'papr_ccdf.png'}")


if __name__ == "__main__":
    main()
