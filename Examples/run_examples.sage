r"""
Collect every example extracted from the notebooks, and time or profile them.

Run from this directory, as the library and the ``tests_*.sage`` files are
loaded by relative path::

    cd Examples

    # list what is available
    sage run_examples.sage --list

    # time a selection (shell glob syntax on the example names)
    sage run_examples.sage 'chen/*' 'merikoski/F[123]'

    # profile one example, in a single process so that cProfile sees what
    # the worker pool would otherwise hide
    sage run_examples.sage --profile --limit 30 maynard/I1

From a Sage session started here::

    sage: load("run_examples.sage")
    sage: run("stadlmann/1*")
    sage: profile("maynard/I1", precision=25)

Each example is a callable taking ``precision`` and ``verbose`` keyword
arguments; the default precision is the one its source notebook uses. A bare
run leaves out ``SLOW_EXAMPLES``, which have to be named explicitly.

Note on parallelism: ``sieve_integral`` and ``sieve_integral_polyratio`` hand
the pieces of a balanced polytope to a ``multiprocessing.Pool`` as soon as
there are more than twenty of them, and cProfile only ever sees the parent
process. ``--profile`` therefore swaps the pool for an in-process map, unless
``--parallel`` is given; wall-clock timings taken the two ways are not
comparable.
"""

EXAMPLES_AUTORUN = False

load("tests_chen.sage")
load("tests_ford_maynard.sage")
load("tests_maynard.sage")
load("tests_merikoski.sage")
load("tests_stadlmann.sage")

import cProfile
import fnmatch
import os
import pstats
import sys
from contextlib import contextmanager
from time import time


ALL_EXAMPLES = {}
for _registry in (CHEN_EXAMPLES, FORD_MAYNARD_EXAMPLES, MAYNARD_EXAMPLES,
                  MERIKOSKI_EXAMPLES, STADLMANN_EXAMPLES):
    ALL_EXAMPLES.update(_registry)


# Examples that are much slower than the rest: excluded from a bare run over
# everything, and reachable by naming them explicitly.
SLOW_EXAMPLES = {f"stadlmann/{j:02d}" for j in STADLMANN_SLOW}


class SerialPool:
    """
    A stand-in for ``multiprocessing.Pool`` that maps in this process.

    The library uses the pool as ``with pool(n) as p: p.map(f, xs)`` and
    nothing else, so that is all the surface this has to reproduce.
    """
    def __init__(self, processes=None):
        pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, func, iterable):
        return [func(x) for x in iterable]


@contextmanager
def serial_workers():
    """
    Run the body with the worker pool replaced by an in-process map.

    The library is brought in with ``load``, so its functions look ``pool`` up
    in this very namespace; rebinding the name here is enough to redirect
    them.
    """
    saved = globals()["pool"]
    globals()["pool"] = SerialPool
    try:
        yield
    finally:
        globals()["pool"] = saved


def select(patterns=("*",)):
    """
    Return the names of the examples matching any of ``patterns``, in
    registration order. A bare ``"*"`` leaves out ``SLOW_EXAMPLES``.
    """
    if isinstance(patterns, str):
        patterns = (patterns,)
    patterns = tuple(patterns)
    names = [name for name in ALL_EXAMPLES
             if any(fnmatch.fnmatch(name, p) for p in patterns)]
    if patterns == ("*",):
        names = [name for name in names if name not in SLOW_EXAMPLES]
    return names


def run(patterns=("*",), precision=None, verbose=0, repeat=1, serial=False):
    """
    Time every matching example and return ``{name: (value, [time, ...])}``.

    ``precision``, if given, overrides each example's own default. An example
    that raises is reported and does not stop the batch; its entry holds the
    exception instead of a value.
    """
    kwds = {"verbose": verbose}
    if precision is not None:
        kwds["precision"] = precision

    results = {}
    for name in select(patterns):
        test = ALL_EXAMPLES[name]
        times = []
        value = None
        try:
            for _ in range(repeat):
                start = time()
                if serial:
                    with serial_workers():
                        value = test(**kwds)
                else:
                    value = test(**kwds)
                times.append(time() - start)
        except Exception as exc:
            print(f"{name:24s} {'--':>8s}   FAILED: "
                  f"{type(exc).__name__}: {exc}")
            results[name] = (exc, times)
            continue
        results[name] = (value, times)
        print(f"{name:24s} {min(times):8.2f}s  {value}")
    total = sum(min(times) for _, times in results.values() if times)
    print(f"{'total (best of each)':24s} {total:8.2f}s")
    return results


def profile(patterns, precision=None, verbose=0, serial=True,
            sort="cumulative", limit=40, dump=None):
    """
    Profile the matching examples under cProfile, print the top ``limit``
    lines of the statistics, and return the ``pstats.Stats``.

    ``serial`` defaults to True here: with the worker pool in place the
    profile would only account for the parent process.
    """
    kwds = {"verbose": verbose}
    if precision is not None:
        kwds["precision"] = precision

    names = select(patterns)
    profiler = cProfile.Profile()

    def body():
        for name in names:
            value = ALL_EXAMPLES[name](**kwds)
            print(f"{name:24s} {value}")

    start = time()
    if serial:
        with serial_workers():
            profiler.runcall(body)
    else:
        profiler.runcall(body)
    elapsed = time() - start

    print(f"\nprofiled {len(names)} example(s) in {elapsed:.2f}s"
          f"{' (serial)' if serial else ' (parallel)'}\n")

    stats = pstats.Stats(profiler)
    if dump is not None:
        stats.dump_stats(dump)
        print(f"raw profile written to {dump}\n")
    stats.sort_stats(sort).print_stats(limit)
    return stats


def main(argv):
    import argparse

    parser = argparse.ArgumentParser(
        prog="run_examples.sage",
        description="Time or profile the sieve integral examples extracted "
                    "from the notebooks.")
    parser.add_argument("patterns", nargs="*", default=["*"],
                        help="shell globs on example names "
                             "(default: everything but the slow ones)")
    parser.add_argument("--list", action="store_true",
                        help="list the matching example names and exit")
    parser.add_argument("--profile", action="store_true",
                        help="run under cProfile instead of just timing")
    parser.add_argument("--serial", action="store_true",
                        help="replace the worker pool by an in-process map "
                             "(implied by --profile)")
    parser.add_argument("--parallel", action="store_true",
                        help="keep the worker pool even under --profile, "
                             "which then only sees the parent process")
    parser.add_argument("--precision", type=int, default=None,
                        help="override each example's default precision")
    parser.add_argument("--repeat", type=int, default=1,
                        help="number of timing repetitions (default 1)")
    parser.add_argument("--verbose", type=int, default=0,
                        help="verbosity passed to the library")
    parser.add_argument("--sort", default="cumulative",
                        help="pstats sort key (default: cumulative)")
    parser.add_argument("--limit", type=int, default=40,
                        help="number of profile lines to print")
    parser.add_argument("--dump", default=None,
                        help="write the raw profile to this file")
    args = parser.parse_args(argv)

    patterns = tuple(args.patterns) or ("*",)

    if args.list:
        for name in select(patterns):
            print(name)
        return

    if args.profile:
        profile(patterns,
                precision=args.precision,
                verbose=args.verbose,
                serial=not args.parallel,
                sort=args.sort,
                limit=args.limit,
                dump=args.dump)
    else:
        run(patterns,
            precision=args.precision,
            verbose=args.verbose,
            repeat=args.repeat,
            serial=args.serial)


# ``__name__`` is "__main__" both when this file is run as a script and when it
# is ``load``ed into a Sage session, so it cannot tell the two apart: guarding
# on it alone made ``load("run_examples.sage")``, which this file's own
# docstring suggests, compute all eighty examples before returning.
# ``sys.argv[0]`` names the preparsed script only in the first case.
if os.path.basename(sys.argv[0]).startswith("run_examples"):
    main(sys.argv[1:])
