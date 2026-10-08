"""Command line: `uv run aqf <step>`. Each step writes to outputs/analysis/<NN_step>/."""

from __future__ import annotations

import argparse
import logging

from aqf.data import build_s5p_cube, load_config

STEPS = ["inventory", "cube", "quality", "trends", "spatial", "s2"]


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="aqf", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("inventory", help="manifest/file integrity, grids, SHA-256 provenance, availability")
    p.add_argument("--no-hash", action="store_true", help="skip SHA-256 (faster)")
    p = sub.add_parser("cube", help="build data/processed/s5p_cube.npz")
    p.add_argument("--force", action="store_true")
    sub.add_parser("quality", help="S5P coverage, zero-floor and outlier analysis")
    sub.add_parser("trends", help="time series, seasonality, STL, Mann-Kendall trends, lockdown")
    sub.add_parser("spatial", help="mean/seasonal maps, hotspots, industrial hubs")
    sub.add_parser("s2", help="Sentinel-2 coverage and spectral-index summary (streamed)")
    p = sub.add_parser("analysis", help="run every step above in order")
    p.add_argument("--no-hash", action="store_true")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    # GDAL warns about a photometric/ExtraSamples TIFF header quirk in the S2 files; band data is unaffected
    logging.getLogger("rasterio").setLevel(logging.ERROR)

    cfg = load_config()
    run_all = args.cmd == "analysis"
    if args.cmd == "inventory" or run_all:
        from aqf.analysis import inventory

        inventory.run(cfg, hash_files=not args.no_hash)
    if args.cmd in ("cube", "quality", "trends", "spatial") or run_all:
        cube = build_s5p_cube(cfg, force=getattr(args, "force", False) or run_all)
    if args.cmd == "quality" or run_all:
        from aqf.analysis import quality

        quality.run(cube, cfg)
    if args.cmd == "trends" or run_all:
        from aqf.analysis import trends

        trends.run(cube, cfg)
    if args.cmd == "spatial" or run_all:
        from aqf.analysis import spatial

        spatial.run(cube, cfg)
    if args.cmd == "s2" or run_all:
        from aqf.analysis import s2

        s2.run(cfg)
