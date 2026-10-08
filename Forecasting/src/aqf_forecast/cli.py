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
    sub.add_parser("ablations", help="F4: ablations, rolling-origin years 2021-2023, NE-monsoon gate (needs F2 first)")
    q = sub.add_parser("final", help="F5: score the 2024 test ONCE with the frozen choice, then the Jan-May 2025 outlook")
    q.add_argument("--rescore-test", action="store_true",
                   help="re-run the 2024 scoring (only to reproduce the identical numbers; never to tune)")
    sub.add_parser("outlook", help="F5 outlook only: refit on 2019-2024 and forecast Jan-May 2025")
    sub.add_parser("all", help="F0, F1, F2 and F4 (F5 is run separately: `final`)")
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
    if args.cmd in ("ablations", "all"):
        from aqf_forecast import ablations

        ablations.run(p)
    if args.cmd in ("final", "outlook"):
        from aqf_forecast import final

        if args.cmd == "final":
            final.test(p, rescore=args.rescore_test)
        final.outlook(p)
