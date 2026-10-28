"""Command line interface: `cp2-analysis <command> [inputs] [options]`."""

import argparse

from .commands import acf, check_exact, raw, summary, thermalisation

COMMANDS = {
    "summary": (summary, "magnetisation and susceptibility vs the exact results"),
    "raw": (raw, "block-averaged raw columns (thermalisation / autocorrelation by eye)"),
    "acf": (acf, "autocorrelation function and integrated autocorrelation time"),
    "thermalisation": (thermalisation, "skip scan and MSER: how many measurements to discard"),
    "check-exact": (check_exact, "pass/fail: M and chi within --max-pull errors of the exact results"),
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cp2-analysis", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    for name, (module, help_text) in COMMANDS.items():
        p = sub.add_parser(name, help=help_text, description=module.__doc__)
        module.add_arguments(p)
        p.set_defaults(func=module.run)
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    args.func(args)
