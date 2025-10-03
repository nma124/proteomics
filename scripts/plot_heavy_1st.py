#!/usr/bin/env python3
"""Generate calibration plots for the heavy_1st workflow output."""

from __future__ import annotations

import argparse
import os
import pathlib

_DEFAULT_MPLCONFIG = pathlib.Path('reports/mpl_cache')
if 'MPLCONFIGDIR' not in os.environ:
    _DEFAULT_MPLCONFIG.mkdir(parents=True, exist_ok=True)
    os.environ['MPLCONFIGDIR'] = str(_DEFAULT_MPLCONFIG.resolve())

if 'XDG_CACHE_HOME' not in os.environ:
    xdg_cache = _DEFAULT_MPLCONFIG / 'xdg_cache'
    (xdg_cache / 'fontconfig').mkdir(parents=True, exist_ok=True)
    os.environ['XDG_CACHE_HOME'] = str(xdg_cache.resolve())

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import linear_model

TOL = 1e-6


def _clean_fragment_df(cat_df: pd.DataFrame) -> pd.DataFrame:
    """Drop rows with invalid concentrations/ratios before fitting."""
    cleaned = cat_df.replace([np.inf, -np.inf], np.nan).dropna(subset=["heavy_conc", "area_ratio"])
    if "heavy_conc" in cleaned:
        cleaned = cleaned.loc[~np.isclose(cleaned["heavy_conc"], 0.0, atol=TOL)]
    return cleaned


def _fit_line(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Return fitted y values and model parameters."""
    weights = np.where(np.abs(x) > TOL, 1.0 / x, 1.0)
    model = linear_model.LinearRegression()
    model.fit(x[:, np.newaxis], y[:, np.newaxis], sample_weight=weights)
    y_fit = model.predict(x[:, np.newaxis]).squeeze()
    intercept = float(model.intercept_.squeeze())
    slope = float(model.coef_.squeeze())
    return y_fit, intercept, slope


def plot_fragment(cat_df: pd.DataFrame, title: str, out_path: pathlib.Path) -> None:
    """Produce a scatter + linear-fit plot for a fragment."""
    cleaned = _clean_fragment_df(cat_df)
    if cleaned.empty:
        return

    x = cleaned["heavy_conc"].to_numpy(dtype=float)
    y = cleaned["area_ratio"].to_numpy(dtype=float)
    order = np.argsort(x)
    x_sorted, y_sorted = x[order], y[order]

    y_fit, intercept, slope = _fit_line(x_sorted, y_sorted)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x_sorted, y_sorted, color="C0", label="observed")
    ax.plot(x_sorted, y_fit, color="C1", label="linear fit")
    ax.set_title(title)
    ax.set_xlabel("Heavy concentration (ng/mL)")
    ax.set_ylabel("Area ratio (light / heavy)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    ax.annotate(f"slope={slope:.3g}\nintercept={intercept:.3g}",
                xy=(0.05, 0.95), xycoords="axes fraction", va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def plot_group(cat_df: pd.DataFrame, title: str, out_path: pathlib.Path) -> None:
    """Plot combined fragments belonging to the same peptide group."""
    cleaned = _clean_fragment_df(cat_df)
    if cleaned.empty:
        return

    # Aggregate all fragment data for the group.
    x = cleaned["heavy_conc"].to_numpy(dtype=float)
    y = cleaned["area_ratio"].to_numpy(dtype=float)
    order = np.argsort(x)
    x_sorted, y_sorted = x[order], y[order]

    y_fit, intercept, slope = _fit_line(x_sorted, y_sorted)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(x_sorted, y_sorted, color="C0", label="observed")
    ax.plot(x_sorted, y_fit, color="C1", label="linear fit")
    ax.set_title(title)
    ax.set_xlabel("Heavy concentration (ng/mL)")
    ax.set_ylabel("Area ratio (light / heavy)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.25)
    ax.annotate(f"slope={slope:.3g}\nintercept={intercept:.3g}",
                xy=(0.05, 0.95), xycoords="axes fraction", va="top")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)


def compute_plot_categories(df: pd.DataFrame) -> pd.DataFrame:
    """Ensure plot category columns exist for downstream plotting."""
    if "plot_cat" not in df.columns:
        df["plot_cat"] = df.apply(
            lambda row: f"{row['peptide']}_{row['fragment ion']}", axis=1
        )
    if "plot_cat_grp" not in df.columns:
        df["plot_cat_grp"] = df["plot_cat"].apply(lambda x: "_".join(x.split("_")[:-1]))
    return df


def generate_plots(result_path: pathlib.Path, output_dir: pathlib.Path) -> None:
    """Load heavy_1st results and create per-fragment & grouped plots."""
    df = pd.read_csv(result_path)
    df.columns = df.columns.str.strip().str.lower()

    required_cols = {"peptide", "fragment ion", "heavy_conc", "area_ratio"}
    missing_cols = required_cols - set(df.columns)
    if missing_cols:
        raise ValueError(f"Missing required columns: {sorted(missing_cols)}")

    df = compute_plot_categories(df)

    frag_dir = output_dir / "fragments"
    for plot_cat, cat_df in df.groupby("plot_cat"):
        fragment_name = plot_cat.replace("/", "-")
        out_path = frag_dir / f"{fragment_name}.png"
        plot_fragment(cat_df, plot_cat, out_path)

    grp_dir = output_dir / "groups"
    for grp_cat, grp_df in df.groupby("plot_cat_grp"):
        group_name = grp_cat.replace("/", "-")
        out_path = grp_dir / f"{group_name}.png"
        plot_group(grp_df, grp_cat, out_path)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        type=pathlib.Path,
        default=pathlib.Path("data/output/heavy_1st_workflow_results.csv"),
        help="Path to workflow results CSV.",
    )
    parser.add_argument(
        "--out",
        type=pathlib.Path,
        default=pathlib.Path("reports/heavy_1st/plots"),
        help="Directory to save generated plots.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    generate_plots(args.results, args.out)


if __name__ == "__main__":
    main()
