"""Visualise any GeoTIFF from this project on its own: raw Sentinel-5P, raw Sentinel-2, or the 2025 outlook files.

Usage (from Scratch/):
    uv run python tools/view_tiff.py FILE.tif                    # summary + a grid of every band (up to --max-bands)
    uv run python tools/view_tiff.py FILE.tif --bands 1 3        # only these bands (1-based)
    uv run python tools/view_tiff.py FILE.tif --rgb              # Sentinel-2 true colour (bands B4, B3, B2)
    uv run python tools/view_tiff.py FILE.tif --info             # print metadata only, no figure
    uv run python tools/view_tiff.py FILE.tif --out fig.png      # save instead of opening a window

Examples:
    uv run python tools/view_tiff.py Dataset/data/s5p_composites/s5p_20200301.tif
    uv run python tools/view_tiff.py Dataset/data/s2_composites/s2_20200301.tif --rgb
    uv run python tools/view_tiff.py Forecasting/results/F5_final/outlook_2025/outlook_NO2_2025.tif --bands 1 6 12 18 24 30

Nodata (-inf / NaN) is shown in grey. Colours are stretched between the 2nd and 98th percentile of each band.
"""

from __future__ import annotations

import argparse
import logging
import math
from pathlib import Path

import numpy as np
import rasterio

S2_NAMES = ["B2 blue", "B3 green", "B4 red", "B5 red-edge 1", "B6 red-edge 2", "B7 red-edge 3", "B8 NIR",
            "B11 SWIR-1", "B12 SWIR-2", "NDVI", "NDBI", "NDMI"]
S5P_NAMES = ["NO2 (mol/m²)", "CO (mol/m²)", "SO2 (mol/m²)"]


def band_names(src) -> list[str]:
    """Band labels: the file's own descriptions if present, else the known S5P / S2 layout, else 'band N'."""
    desc = list(src.descriptions)
    if any(desc):
        return [d or f"band {i + 1}" for i, d in enumerate(desc)]
    if src.count == 3 and Path(src.name).name.startswith("s5p"):
        return S5P_NAMES
    if src.count == 12 and Path(src.name).name.startswith("s2"):
        return S2_NAMES
    return [f"band {i + 1}" for i in range(src.count)]


def read(src, bands: list[int]) -> np.ndarray:
    a = src.read(bands).astype("float64")
    a[~np.isfinite(a)] = np.nan
    if src.nodata is not None and np.isfinite(src.nodata):
        a[a == src.nodata] = np.nan
    return a


def info(src) -> None:
    names = band_names(src)
    b = src.bounds
    print(f"file      : {src.name}")
    print(f"size      : {src.count} band(s) x {src.height} rows x {src.width} cols, {src.dtypes[0]}")
    print(f"CRS       : {src.crs}   pixel size: {src.res[0]:.5f} x {src.res[1]:.5f}")
    print(f"bounds    : lon {b.left:.3f}..{b.right:.3f}, lat {b.bottom:.3f}..{b.top:.3f}")
    if src.tags():
        print(f"file tags : {src.tags()}")
    print("bands     :")
    for i in range(1, src.count + 1):
        a = read(src, [i])[0]
        ok = np.isfinite(a)
        stats = (f"min {np.nanmin(a):.4g}  mean {np.nanmean(a):.4g}  max {np.nanmax(a):.4g}" if ok.any()
                 else "all nodata")
        tags = src.tags(i)
        extra = f"  [model: {tags['model']}]" if "model" in tags else ""
        print(f"  {i:>3}  {names[i - 1]:<40} valid {ok.mean():6.1%}  {stats}{extra}")


def plot(src, bands: list[int], rgb: bool, out: str | None) -> None:
    import matplotlib

    if out:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    b = src.bounds
    extent = [b.left, b.right, b.bottom, b.top]
    names = band_names(src)
    if rgb:
        if src.count < 3:
            raise SystemExit("--rgb needs a multi-band image (a Sentinel-2 file)")
        idx = [3, 2, 1] if src.count == 12 else [1, 2, 3]          # S2: B4 red, B3 green, B2 blue
        a = read(src, idx)
        lo, hi = np.nanpercentile(a, [2, 98])
        img = np.clip((a - lo) / (hi - lo), 0, 1).transpose(1, 2, 0)
        img = np.where(np.isfinite(img), img, 0.75)                # nodata in grey
        fig, ax = plt.subplots(figsize=(7, 7))
        ax.imshow(img, extent=extent)
        ax.set_title(f"{Path(src.name).name}: true colour (bands {idx})")
        ax.set_xlabel("lon (°E)")
        ax.set_ylabel("lat (°N)")
    else:
        n = len(bands)
        cols = min(n, 4)
        rows = math.ceil(n / cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4.2 * cols, 3.9 * rows), squeeze=False)
        cmap = plt.get_cmap("viridis").copy()
        cmap.set_bad("#bdbdbd")
        for k, ax in enumerate(axes.ravel()):
            if k >= n:
                ax.axis("off")
                continue
            a = read(src, [bands[k]])[0]
            if np.isfinite(a).any():
                vmin, vmax = np.nanpercentile(a, [2, 98])
            else:
                vmin = vmax = None
            im = ax.imshow(np.ma.masked_invalid(a), extent=extent, cmap=cmap, vmin=vmin, vmax=vmax,
                           interpolation="nearest")
            ax.set_title(f"{bands[k]}: {names[bands[k] - 1]}", fontsize=8)
            ax.tick_params(labelsize=7)
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03).ax.tick_params(labelsize=7)
        fig.suptitle(Path(src.name).name)
        fig.tight_layout()
    if out:
        fig.savefig(out, dpi=150, bbox_inches="tight")
        print(f"saved {out}")
    else:
        plt.show()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("file", help="path to a .tif / .tiff file")
    ap.add_argument("--bands", type=int, nargs="+", help="1-based band numbers to plot (default: all, up to --max-bands)")
    ap.add_argument("--max-bands", type=int, default=12, help="cap on bands plotted when --bands is not given")
    ap.add_argument("--rgb", action="store_true", help="true-colour composite (Sentinel-2)")
    ap.add_argument("--info", action="store_true", help="print metadata and band statistics only")
    ap.add_argument("--out", help="save the figure to this PNG instead of opening a window")
    args = ap.parse_args()
    # GDAL warns about a harmless photometric/ExtraSamples header quirk in the Sentinel-2 files
    logging.getLogger("rasterio").setLevel(logging.ERROR)
    with rasterio.open(args.file) as src:
        info(src)
        if args.info:
            return
        bands = args.bands or list(range(1, min(src.count, args.max_bands) + 1))
        bad = [x for x in bands if not 1 <= x <= src.count]
        if bad:
            raise SystemExit(f"band(s) {bad} out of range 1..{src.count}")
        if not args.bands and src.count > args.max_bands:
            print(f"(showing bands 1-{args.max_bands} of {src.count}; use --bands to choose others)")
        plot(src, bands, args.rgb, args.out)


if __name__ == "__main__":
    main()
