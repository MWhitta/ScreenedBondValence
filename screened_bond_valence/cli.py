"""Command-line entrypoints for ScreenedBondValence."""

from __future__ import annotations

import argparse
from collections.abc import Sequence

from .constants import DEFAULT_ALGORITHMS


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="screened-bond-valence")
    subparsers = parser.add_subparsers(dest="command", required=True)

    batch_parser = subparsers.add_parser("batch", help="Fetch Materials Project data and fit one or more cation-anion systems")
    batch_parser.add_argument("--api-key", required=True, help="Materials Project API key")
    batch_parser.add_argument("--cations", nargs="+", required=True, help="One or more cation symbols")
    batch_parser.add_argument("--anions", nargs="+", required=True, help="One or more anion symbols")
    batch_parser.add_argument(
        "--algorithms",
        nargs="+",
        default=list(DEFAULT_ALGORITHMS),
        help="Optimization algorithms to run",
    )
    batch_parser.add_argument("--output-dir", default="res", help="Directory for fitted outputs")
    batch_parser.add_argument(
        "--r0-bounds",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(0.0, 5.0),
        help="Bounds for R0 optimization",
    )
    batch_parser.add_argument(
        "--energy-above-hull",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(0.0, 0.05),
        help="Materials Project energy-above-hull filter",
    )

    web_parser = subparsers.add_parser("web", help="Launch the Gradio CIF fitting app")
    web_parser.add_argument("--server-name", default="127.0.0.1", help="Bind address for the web app")
    web_parser.add_argument("--server-port", type=int, default=7860, help="Port for the web app")
    web_parser.add_argument("--share", action="store_true", help="Enable Gradio share links")
    web_parser.add_argument("--inbrowser", action="store_true", help="Open the app in a browser")

    return parser


def _build_processor(**processor_kwargs: object):
    from .processor import BondValenceProcessor

    return BondValenceProcessor(**processor_kwargs)


def run_batch(args: argparse.Namespace) -> int:
    processor = _build_processor(
        api_key=args.api_key,
        algos=args.algorithms,
        cations=args.cations,
        anions=args.anions,
        output_dir=args.output_dir,
        r0_bounds=tuple(args.r0_bounds),
        energy_above_hull=tuple(args.energy_above_hull),
    )
    for cation in args.cations:
        for anion in args.anions:
            processor.process_cation_system(cation, anion)
    return 0


def _launch_web_app(**launch_kwargs: object) -> None:
    from .web import launch_app

    launch_app(**launch_kwargs)


def run_web(args: argparse.Namespace) -> int:
    _launch_web_app(
        server_name=args.server_name,
        server_port=args.server_port,
        share=args.share,
        inbrowser=args.inbrowser,
    )
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(list(argv) if argv is not None else None)

    if args.command == "batch":
        return run_batch(args)
    if args.command == "web":
        return run_web(args)

    parser.error(f"Unknown command: {args.command}")
    return 2
