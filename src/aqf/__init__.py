"""aqf: satellite air-pollution analysis and forecasting over Chennai (Sentinel-5P + Sentinel-2)."""


def main() -> None:
    from aqf.cli import main as _main

    _main()
