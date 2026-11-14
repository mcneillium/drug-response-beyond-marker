"""Figures for the study, drawn from the tables in results/.

    python scripts/plot_results.py

    figures/ladder.png                  what each kind of information adds
    figures/leak.png                    feature selection done wrongly and correctly
    figures/overlap.png                 four ways to handle missing protein data
    figures/cluster_stability.png       cluster stability against a no-groups reference
    figures/explanation_stability.png   how often each feature is in a model's top 20
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SURFACE, TEXT, TEXT_QUIET, GRID = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df"
BLUE, ORANGE, GREEN, YELLOW, NEUTRAL = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#a8a7a2"


def tidy_axis(axis, grid: str = "y") -> None:
    axis.set_facecolor(SURFACE)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    axis.spines["bottom"].set_color(GRID)
    axis.tick_params(colors=TEXT_QUIET, labelsize=8.5, length=0)
    (axis.yaxis if grid == "y" else axis.xaxis).grid(True, color=GRID, linewidth=0.8)
    axis.set_axisbelow(True)


def save(figure, path: Path) -> None:
    figure.savefig(path, dpi=200, facecolor=SURFACE, bbox_inches="tight")
    plt.close(figure)


def ladder_figure(table: pd.DataFrame, path: Path) -> None:
    steps = ["L1", "L2", "L3", "L4", "L5"]
    names = ["PTEN status\nalone", "Cancer type +\nsequencing centre", "PTEN +\nbackground",
             "+ gene expression\n(elastic net)", "+ gene expression\n(boosting)"]
    rows = table.set_index("step").loc[steps]
    figure, axis = plt.subplots(figsize=(8.2, 4.2), facecolor=SURFACE)
    bars = axis.bar(range(len(steps)), rows["r2_mean"], width=0.55, color=BLUE, edgecolor=SURFACE, linewidth=2)
    axis.errorbar(range(len(steps)), rows["r2_mean"], yerr=rows["r2_sd"], fmt="none", ecolor=TEXT, elinewidth=1.2,
                  capsize=3)
    for bar, (_, row) in zip(bars, rows.iterrows()):
        axis.annotate(f"{row['r2_mean']:.2f}", (bar.get_x() + bar.get_width() / 2, row["r2_mean"] + row["r2_sd"]),
                      xytext=(0, 5), textcoords="offset points", ha="center", fontsize=9, color=TEXT)
    tidy_axis(axis)
    axis.set_xticks(range(len(steps)), names)
    axis.set_ylim(0, 0.30)
    axis.set_ylabel("Share of the variation in response explained (R²)", fontsize=8.5, color=TEXT_QUIET)
    axis.set_title("More than four fifths of the response stays unexplained at every step",
                   fontsize=11.5, color=TEXT, loc="left", pad=22)
    axis.text(0, 1.035, "Ipatasertib, 924 cell lines, held-out predictions. Lines show the spread over 5 repeats.",
              transform=axis.transAxes, fontsize=8.5, color=TEXT_QUIET)
    save(figure, path)


def leak_figure(table: pd.DataFrame, path: Path) -> None:
    order = [("real response", "wrong"), ("real response", "right"), ("shuffled response", "wrong"),
             ("shuffled response", "right")]
    values, spreads, colours = [], [], []
    for response, kind in order:
        row = table[(table["response"] == response) & table["method"].str.startswith(kind)].iloc[0]
        values.append(row["r2_mean"]), spreads.append(row["r2_sd"])
        colours.append(ORANGE if kind == "wrong" else BLUE)
    positions = [0, 1, 2.6, 3.6]
    figure, axis = plt.subplots(figsize=(7.6, 4.2), facecolor=SURFACE)
    axis.bar(positions, values, width=0.7, color=colours, edgecolor=SURFACE, linewidth=2)
    axis.errorbar(positions, values, yerr=spreads, fmt="none", ecolor=TEXT, elinewidth=1.2, capsize=3)
    for x, value, spread in zip(positions, values, spreads):
        shown = round(value, 2) + 0.0                    # avoids printing "-0.00"
        above = shown >= 0
        axis.annotate(f"{shown:+.2f}" if shown else "0.00", (x, value + spread if above else value - spread),
                      xytext=(0, 5) if above else (0, -13), textcoords="offset points", ha="center",
                      fontsize=9, color=TEXT)
    axis.axhline(0, color=TEXT_QUIET, linewidth=0.8)
    tidy_axis(axis)
    axis.set_xticks([0.5, 3.1], ["Real drug response", "Shuffled response\n(nothing real to find)"])
    axis.tick_params(axis="x", labelsize=9.5)
    axis.set_ylim(-0.27, 0.30)
    axis.set_ylabel("Held-out R²", fontsize=8.5, color=TEXT_QUIET)
    handles = [plt.Rectangle((0, 0), 1, 1, color=ORANGE), plt.Rectangle((0, 0), 1, 1, color=BLUE)]
    axis.legend(handles, ["Wrong: genes chosen using all cell lines", "Right: genes chosen inside each training fold"],
                loc="upper right", fontsize=8.5, frameon=False, labelcolor=TEXT)
    axis.set_title("Choosing genes before cross-validating invents skill that is not there",
                   fontsize=11.5, color=TEXT, loc="left", pad=22)
    axis.text(0, 1.035, "200 genes most correlated with response, ridge regression. Lines show the spread over 5 repeats.",
              transform=axis.transAxes, fontsize=8.5, color=TEXT_QUIET)
    save(figure, path)


def overlap_figure(table: pd.DataFrame, path: Path) -> None:
    series = {"expression_only": ("Expression only", BLUE, "o"), "late_fusion": ("Late fusion", GREEN, "^"),
              "fill_gaps": ("Fill the gaps", ORANGE, "s"), "complete_cases": ("Complete cases", YELLOW, "D")}
    shares = [0.97, 0.75, 0.50, 0.25]
    panels = [("r2_all", "All 924 cell lines", ["expression_only", "late_fusion", "fill_gaps"]),
              ("r2_with_protein", "Only the cell lines that have protein data", list(series))]
    figure, axes = plt.subplots(1, 2, figsize=(10.2, 4.3), facecolor=SURFACE, sharey=True)
    for axis, (column, title, names) in zip(axes, panels):
        for name in names:
            label, colour, marker = series[name]
            rows = table[table["strategy"] == name].set_index("share_with_protein").reindex(shares)
            axis.plot(range(len(shares)), rows[f"{column}_mean"], color=colour, linewidth=2, marker=marker,
                      markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.2, label=label)
        tidy_axis(axis)
        axis.set_xticks(range(len(shares)), [f"{int(round(100 * s))}%" for s in shares])
        axis.set_xlabel("Share of cell lines that keep their protein data", fontsize=8.5, color=TEXT_QUIET)
        axis.set_ylim(0, 0.21)
        axis.set_yticks([0, 0.05, 0.10, 0.15, 0.20])
        axis.set_title(title, fontsize=10.5, color=TEXT, loc="left")
    axes[0].set_ylabel("Held-out R², mean of 5 repeats", fontsize=8.5, color=TEXT_QUIET)
    axes[1].legend(loc="lower left", fontsize=8.5, frameon=False, labelcolor=TEXT)
    figure.suptitle("Protein data did not help here, and two common ways of using it made things worse",
                    fontsize=12, color=TEXT, x=0.065, ha="left", y=0.99)
    save(figure, path)


def cluster_figure(table: pd.DataFrame, path: Path) -> None:
    versions = [("as_published", "Expression as published"), ("adjusted", "Tissue and sequencing centre removed")]
    figure, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), facecolor=SURFACE, sharey=True)
    for axis, (version, title) in zip(axes, versions):
        rows = table[table["version"] == version].sort_values("k")
        ks = rows["k"].to_numpy()
        axis.fill_between(ks, 0, rows["reference_max_ari"], color=NEUTRAL, alpha=0.30, linewidth=0,
                          label="Reference: data with no groups (up to the highest of 20)")
        axis.plot(ks, rows["reference_mean_ari"], color=TEXT_QUIET, linewidth=1.5, linestyle=(0, (4, 3)),
                  label="Reference: mean")
        axis.plot(ks, rows["stability_mean_ari"], color=BLUE, linewidth=2, marker="o", markersize=7,
                  markeredgecolor=SURFACE, markeredgewidth=1.2, label="PTEN-altered cell lines")
        for k, value in zip(ks, rows["stability_mean_ari"]):
            axis.annotate(f"{value:.2f}", (k, value), xytext=(0, 10), textcoords="offset points", ha="center",
                          fontsize=8.5, color=TEXT,
                          bbox={"facecolor": SURFACE, "edgecolor": "none", "pad": 1.2, "alpha": 0.85})
        tidy_axis(axis)
        axis.set_xticks(ks, [f"{k} groups" for k in ks])
        axis.set_ylim(0, 1.12)
        axis.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        axis.set_title(title, fontsize=10.5, color=TEXT, loc="left")
    axes[0].set_ylabel("Stability (mean adjusted Rand index)", fontsize=8.5, color=TEXT_QUIET)
    axes[1].legend(loc="lower left", fontsize=8, frameon=False, labelcolor=TEXT)
    figure.suptitle("Stable groups exist only before tissue is removed, and they follow tissue",
                    fontsize=12, color=TEXT, x=0.065, ha="left", y=0.99)
    save(figure, path)


def explanation_figure(stability: pd.DataFrame, path: Path, n_fits: int = 25) -> None:
    rows = stability.sort_values(["times_in_top20", "mean_abs_shap"], ascending=False).head(20).iloc[::-1]
    names = rows["feature"].replace({"PTEN_status": "PTEN status (the marker)"})
    figure, axis = plt.subplots(figsize=(7.4, 5.6), facecolor=SURFACE)
    axis.barh(range(len(rows)), rows["times_in_top20"], height=0.62, color=BLUE, edgecolor=SURFACE, linewidth=1.5)
    for y, value in enumerate(rows["times_in_top20"]):
        axis.annotate(f"{int(value)}", (value, y), xytext=(5, 0), textcoords="offset points", va="center",
                      fontsize=8.5, color=TEXT)
    axis.axvline(0.8 * n_fits, color=TEXT_QUIET, linewidth=1, linestyle=(0, (4, 3)))
    axis.annotate("in 80% of fits", (0.8 * n_fits, len(rows) - 0.4), xytext=(-4, 0), textcoords="offset points",
                  ha="right", va="bottom", fontsize=8.5, color=TEXT_QUIET)
    tidy_axis(axis, grid="x")
    axis.set_yticks(range(len(rows)), names)
    axis.tick_params(axis="y", labelsize=9, colors=TEXT)
    axis.set_xlim(0, n_fits + 2)
    axis.set_ylim(-0.6, len(rows) + 0.4)
    axis.set_xlabel(f"Number of the {n_fits} model fits in which the feature is among the 20 most important",
                    fontsize=8.5, color=TEXT_QUIET)
    axis.set_title("Only two features are reliably in the top 20", fontsize=11.5, color=TEXT, loc="left", pad=22)
    axis.text(0, 1.03, "Gradient boosting, importance by mean absolute SHAP value on held-out cell lines.",
              transform=axis.transAxes, fontsize=8.5, color=TEXT_QUIET)
    save(figure, path)


def main() -> None:
    results, figures = ROOT / "results", ROOT / "figures"
    figures.mkdir(exist_ok=True)
    ladder_figure(pd.read_csv(results / "table_ladder.csv"), figures / "ladder.png")
    leak_figure(pd.read_csv(results / "table_leak.csv"), figures / "leak.png")
    overlap_figure(pd.read_csv(results / "table_overlap.csv"), figures / "overlap.png")
    cluster_figure(pd.read_csv(results / "table_clusters.csv"), figures / "cluster_stability.png")
    explanation_figure(pd.read_csv(results / "table_feature_stability.csv"), figures / "explanation_stability.png")
    print("Figures written to", figures)


if __name__ == "__main__":
    main()
