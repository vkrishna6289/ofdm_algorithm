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
    written = ["ber.png", "papr_ccdf.png"]

    if (out / "pa_psd.csv").exists():
        plot_pa(out)
        written += ["pa_psd.png", "pa_metrics.png", "pa_ber.png"]
    print("Wrote " + ", ".join(str(out / w) for w in written))


def plot_pa(out):
    psd = load(out / "pa_psd.csv")
    k = psd["subcarrier"]
    fig, ax = plt.subplots(figsize=(8, 5))
    for col, style in [("linear_rect", ":"), ("linear", "-"), ("ibo10", "-"), ("ibo6", "-"),
                       ("ibo2", "-"), ("ibo0", "-")]:
        name = {"linear_rect": "linear PA, no TX window", "linear": "linear PA"}.get(col, col.replace("ibo", "IBO ") + " dB")
        ax.plot(k, psd[col], style, label=name, linewidth=1.2)
    for edge in (-90, -38, 38, 90):
        ax.axvline(edge, color="grey", alpha=0.3, linewidth=0.8)
    ax.set(xlabel="Frequency [subcarrier index]", ylabel="PSD rel. in-band [dB]",
           title="Spectral regrowth from Rapp PA (p=2); grey lines = adjacent channels", ylim=(-70, 5))
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "pa_psd.png", dpi=150)

    with open(out / "pa_metrics.csv", newline="") as f:
        rows = [r for r in csv.DictReader(f) if float(r["ibo_db"]) >= 0]
    ibo = [float(r["ibo_db"]) for r in rows]
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.plot(ibo, [float(r["aclr_db"]) for r in rows], "o-", color="C0", label="ACLR")
    ax1.set_xlabel("Input back-off [dB]")
    ax1.set_ylabel("ACLR [dB]", color="C0")
    ax2 = ax1.twinx()
    ax2.plot(ibo, [float(r["evm_pct"]) for r in rows], "s-", color="C3", label="EVM")
    ax2.set_ylabel("EVM [%]", color="C3")
    ax1.set_title("PA distortion vs back-off")
    ax1.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "pa_metrics.png", dpi=150)

    ber = load(out / "pa_ber.csv")
    x = ber["ebn0_db"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.semilogy(*masked(x, ber["theory"]), "k-", label="theory (linear)")
    for col in [c for c in ber if c not in ("ebn0_db", "theory")]:
        ax.semilogy(*masked(x, ber[col]), "o-", markersize=4, label=col.replace("ibo", "IBO ") + (" dB" if col != "linear" else ""))
    ax.set(xlabel="Eb/N0 [dB]", ylabel="BER", title="QPSK-OFDM BER through Rapp PA (AWGN)", ylim=(1e-7, 0.2))
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "pa_ber.png", dpi=150)


if __name__ == "__main__":
    main()
