#!/usr/bin/env python
"""Figure 3 — potential energy during the two minimisations, from the .edr files GROMACS wrote.

Run with the package environment's Python (it has matplotlib) after extracting the curves:

    printf 'Potential\n0\n' | gmx energy -f em_hbonds.edr -o em_vacuum.xvg
    printf 'Potential\n0\n' | gmx energy -f em.edr        -o em_solvated.xvg
    python em_energy.py em_vacuum.xvg em_solvated.xvg energy.jpg

Inputs are copies of <working_dir>/pdbmdauto-e2e-full/e2e_4z8j/gmx/{em_hbonds.edr,em.edr}.
"""
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_xvg(path):
    xs, ys = [], []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(("#", "@")) or not line.strip():
                continue
            parts = line.split()
            xs.append(float(parts[0]))
            ys.append(float(parts[1]))
    return xs, ys


def main(vacuum, solvated, out):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2), dpi=120)
    for ax, (path, title) in zip(
        axes,
        [(vacuum, "In vacuum, restrained (step 6)"), (solvated, "In water, after solvation (step 10)")],
    ):
        xs, ys = read_xvg(path)
        ax.plot(xs, [y / 1000 for y in ys], color="#1f6f78", linewidth=1.6)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel("minimisation step")
        ax.set_ylabel("potential energy (10³ kJ/mol)")
        ax.grid(alpha=0.3)
        ax.annotate(
            f"{ys[-1]/1000:,.1f}",
            xy=(xs[-1], ys[-1] / 1000),
            xytext=(-40, 12),
            textcoords="offset points",
            fontsize=9,
            color="#1f6f78",
        )
    fig.suptitle("Steepest-descent minimisation converges — the number to learn to read", fontsize=12)
    fig.tight_layout()
    fig.savefig(out, format="jpg", pil_kwargs={"quality": 88})
    print(f"wrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
