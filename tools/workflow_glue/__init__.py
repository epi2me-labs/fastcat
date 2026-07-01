"""Vendored wf-template workflow glue used for report construction."""

import argparse
import glob
import importlib
import logging
import os
import sys

from .util import _log_level, get_main_logger


__version__ = "0.0.1"
_package_name = "workflow_glue"


def get_components(allowed_components=None):
    """Find workflow glue commands with argparser() and main()."""
    logger = get_main_logger(_package_name)
    home_path = os.path.dirname(os.path.abspath(__file__))

    components = {}
    for fname in glob.glob(os.path.join(home_path, "*.py")):
        name = os.path.splitext(os.path.basename(fname))[0]
        if name in ("__init__", "util"):
            continue
        if allowed_components is not None and name not in allowed_components:
            continue

        try:
            mod = importlib.import_module(f"{_package_name}.{name}")
        except ModuleNotFoundError as e:
            logger.warning("could not load %s due to missing module %s", name, e.name)
            continue

        if callable(getattr(mod, "main", None)) and callable(getattr(mod, "argparser", None)):
            components[name] = mod
    return components


def cli():
    """Run workflow glue entrypoints."""
    logger = get_main_logger(_package_name)
    parser = argparse.ArgumentParser(
        "dev-tools",
        parents=[_log_level()],
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s {}".format(__version__))

    subparsers = parser.add_subparsers(
        title="subcommands",
        description="valid commands",
        help="additional help",
        dest="command")
    subparsers.required = True

    if len(sys.argv) > 1:
        components = get_components(allowed_components=[sys.argv[1]])
        if sys.argv[1] not in components:
            logger.warning("importing all modules, this may take some time")
            components = get_components()
    else:
        components = get_components()

    for name, module in components.items():
        subparser = subparsers.add_parser(name, parents=[module.argparser()])
        subparser.set_defaults(func=module.main)

    args = parser.parse_args()
    logging.getLogger().setLevel(args.log_level)
    return args.func(args)
