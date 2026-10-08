"""Command line: `uv run aqf-forecast <phase>`. Results go to Forecasting/results/<phase>/."""

from __future__ import annotations

import argparse
import logging


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="aqf-forecast", description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("forecastability", help="F0: autocorrelation, cross-gas, neighbours, S2 vs long-term mean")
    sub.add_parser("baselines", help="F1: climatology / persistence / damped persistence, tasks A and B on validation")
    sub.add_parser("learned", help="F2: LightGBM + Ridge for tasks A and B, compared with the baselines on validation")
    sub.add_parser("all", help="F0, F1 and F2")
    args = ap.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s: %(message)s", datefmt="%H:%M:%S")
    logging.getLogger("rasterio").setLevel(logging.ERROR)   # harmless S2 TIFF header warning

    from aqf_forecast.problem import load_problem

    p = load_problem()
    if args.cmd in ("forecastability", "all"):
        from aqf_forecast import forecastability

        forecastability.run(p)
    if args.cmd in ("baselines", "all"):
        from aqf_forecast import evaluate

        evaluate.run(p, "val")
    if args.cmd in ("learned", "all"):
        from aqf_forecast import learned

        learned.run(p)
