"""Command line: python -m chip_lod {run,figures} [--params P] [--out DIR]."""

import argparse
import time

from .analyses import run_all
from .io import read_params


def main(argv=None):
    ap = argparse.ArgumentParser(prog="chip_lod")
    ap.add_argument("command", choices=["run", "figures"])
    ap.add_argument("--params", default="config/params.json")
    ap.add_argument("--out", default="results/python")
    ap.add_argument("--figdir", default="results/figures")
    args = ap.parse_args(argv)
    if args.command == "run":
        t0 = time.time()
        run_all(read_params(args.params), args.out)
        print(f"wrote results to {args.out} in {time.time() - t0:.1f} s")
    else:
        from .figures import make_all
        make_all(args.out, args.figdir)
        print(f"wrote figures to {args.figdir}")


if __name__ == "__main__":
    main()
