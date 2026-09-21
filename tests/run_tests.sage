r"""
Run the test suite. From the repository root::

    sage tests/run_tests.sage                       # the quick tests
    sage tests/run_tests.sage --long                # quick and long
    sage tests/run_tests.sage --only-timing         # the micro-benchmarks
    sage tests/run_tests.sage 'regression/chen*'    # shell globs on test names

    sage tests/run_tests.sage --long \
        --time "Work plan/runs/04-test-suite-<commit>.json"
    sage tests/run_tests.sage --only-timing \
        --time new.json --compare old.json

Each ``tests/test_*.sage`` file defines

- ``TESTS``: a dict ``name -> (callable, tags)``, the tags taken from
  ``"quick"``, ``"long"`` and ``"timing"``. A test raises ``AssertionError``
  with an informative message when it fails, and returns nothing.
- ``KNOWN_ISSUES`` (optional): a dict ``name -> reason``. Those tests are
  expected to fail; a failure is reported as ``xfail`` and does not make the
  run fail, and a pass is reported as ``XPASS``, which is worth looking at.

Options:

- ``--long`` adds the tests tagged ``long``; ``--only-timing`` runs the
  ``timing`` ones instead, which are otherwise never run;
- ``--time OUT.json`` writes the wall-clock time of every test that ran, with
  the commit, the date and the machine, ``--repeat N`` (default 3, and 1 when
  ``--time`` is absent) taking the median of N runs;
- ``--compare BASE.json`` prints a before/after table against such a file;
- ``--serial`` makes the library map the pieces of a balanced polytope in
  process instead of over a pool of workers. The parallel time is the headline
  number; a serial time is only comparable with another serial time.

Timings mean nothing unless the machine is idle, and a before/after pair
should be taken back to back.
"""
import argparse
import fnmatch
import json
import os
import statistics
import sys
from time import time

load("tests/harness.sage")

TEST_FILES = sorted(f for f in os.listdir("tests")
                    if f.startswith("test_") and f.endswith(".sage"))

ALL_TESTS = {}
ALL_KNOWN_ISSUES = {}
for _file in TEST_FILES:
    TESTS = {}
    KNOWN_ISSUES = {}
    load("tests/" + _file)
    ALL_TESTS.update(TESTS)
    ALL_KNOWN_ISSUES.update(KNOWN_ISSUES)


def select(patterns, tags):
    """
    The tests whose tags meet ``tags``.

    A test tagged ``timing`` is reached only by ``--only-timing``, and one
    tagged ``long`` only when ``--long`` is given, whatever their other tags:
    without the first rule ``--long`` would pull in a benchmark tagged
    ``("timing", "long")``, which is a ten minute one.
    """
    def wanted(test_tags):
        if ("timing" in test_tags) != ("timing" in tags):
            return False
        if "long" in test_tags and "long" not in tags:
            return False
        return bool(set(test_tags) & tags)

    names = [name for name, (_, test_tags) in ALL_TESTS.items()
             if wanted(test_tags)]
    if patterns:
        names = [name for name in names
                 if any(fnmatch.fnmatch(name, p) for p in patterns)]
    return names


def run(names, repeat=1, serial=False, timing=False):
    """
    Run the named tests, and return a list of result dicts.

    When ``timing``, the cached example values are dropped before each run, so
    that every test does its own work and the time recorded is the
    computation's. Without it the cache stands, which is what makes the suite
    quick, and the times printed are not comparable between tests.
    """
    results = []
    width = max(len(name) for name in names) + 2
    for name in names:
        test = ALL_TESTS[name][0]
        times = []
        error = None
        for _ in range(repeat):
            if timing:
                clear_values()
            start = time()
            try:
                if serial:
                    with serial_workers():
                        test()
                else:
                    test()
            except Exception as exc:
                error = f"{type(exc).__name__}: {exc}"
                times.append(time() - start)
                break
            times.append(time() - start)
        # float(): in a .sage file round() is Sage's and returns a Sage real,
        # which json cannot serialise.
        seconds = float(round(statistics.median(times), 3))
        known = ALL_KNOWN_ISSUES.get(name)
        if error is None:
            status = "XPASS" if known else "pass"
        else:
            status = "xfail" if known else "FAIL"
        results.append({"name": name, "status": status,
                        "seconds": seconds, "error": error})
        print(f"{status:>5s}  {name:{width}s} {seconds:7.2f}s", flush=True)
        if status == "FAIL":
            print(f"       {error}", flush=True)
        elif status == "xfail":
            print(f"       known issue: {known}", flush=True)
            print(f"       {error}", flush=True)
        elif status == "XPASS":
            print(f"       known issue, but it passed: {known}", flush=True)
    return results


def summarise(results):
    counts = {}
    for r in results:
        counts[r["status"]] = counts.get(r["status"], 0) + 1
    total = sum(r["seconds"] for r in results)
    print("\n" + "  ".join(f"{n} {s}" for s, n in sorted(counts.items()))
          + f"   in {total:.1f}s")
    return counts


def write_times(path, results, serial, repeat):
    record = {"provenance": provenance(), "serial": bool(serial),
              "repeat": int(repeat),
              "tests": {r["name"]: r["seconds"] for r in results},
              "status": {r["name"]: r["status"] for r in results}}
    with open(path, "w") as f:
        json.dump(record, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"\ntimes written to {path}")


def compare(path, results):
    with open(path) as f:
        base = json.load(f)
    print(f"\nagainst {path} "
          f"({base['provenance']['commit']}, {base['provenance']['date']}, "
          f"{'serial' if base.get('serial') else 'parallel'}):\n")
    print(f"{'test':44s} {'before':>9s} {'after':>9s} {'ratio':>7s}")
    shown = 0
    for r in results:
        before = base["tests"].get(r["name"])
        if before is None:
            continue
        after = r["seconds"]
        ratio = after / before if before else float("nan")
        print(f"{r['name']:44s} {before:8.2f}s {after:8.2f}s {ratio:7.2f}")
        shown += 1
    missing = [r["name"] for r in results if r["name"] not in base["tests"]]
    if missing:
        print(f"\nnot in the baseline: {', '.join(missing)}")
    total_before = sum(base["tests"].get(r["name"], 0) for r in results)
    total_after = sum(r["seconds"] for r in results
                      if r["name"] in base["tests"])
    if shown and total_before:
        print(f"\n{'total of the ' + str(shown) + ' common tests':44s} "
              f"{total_before:8.2f}s {total_after:8.2f}s "
              f"{total_after / total_before:7.2f}")


def main(argv):
    parser = argparse.ArgumentParser(
        prog="run_tests.sage",
        description="Run the sieve integral test suite, from the repository "
                    "root.")
    parser.add_argument("patterns", nargs="*",
                        help="shell globs on test names")
    parser.add_argument("--long", action="store_true",
                        help="also run the tests tagged long")
    parser.add_argument("--only-timing", action="store_true",
                        help="run the micro-benchmarks instead of the tests")
    parser.add_argument("--list", action="store_true",
                        help="list the selected test names and exit")
    parser.add_argument("--time", metavar="OUT.json",
                        help="record the wall-clock times")
    parser.add_argument("--compare", metavar="BASE.json",
                        help="compare the times with an earlier record")
    parser.add_argument("--repeat", type=int, default=None,
                        help="runs per test, median kept (default 3 with "
                             "--time, 1 without)")
    parser.add_argument("--serial", action="store_true",
                        help="map the pieces in process, without the pool")
    args = parser.parse_args(argv)

    tags = {"timing"} if args.only_timing else {"quick"}
    if args.long:
        tags.add("long")
    names = select(args.patterns, tags)

    if args.list:
        for name in names:
            print(name)
        return 0
    if not names:
        print("no test selected")
        return 1

    repeat = args.repeat if args.repeat is not None else (3 if args.time else 1)
    print(f"{len(names)} tests, {'serial' if args.serial else 'parallel'}"
          f"{f', median of {repeat}' if repeat > 1 else ''}\n")
    results = run(names, repeat=repeat, serial=args.serial)
    counts = summarise(results)

    if args.time:
        write_times(args.time, results, args.serial, repeat)
    if args.compare:
        compare(args.compare, results)

    return 1 if counts.get("FAIL") else 0


if os.path.basename(sys.argv[0]).startswith("run_tests"):
    _status = main(sys.argv[1:])
    # sys.exit() does not reach the shell from a .sage script: Sage catches the
    # SystemExit, prints its code and exits 1 whatever it was. os._exit sets the
    # process status, which is what "exits non-zero on failure" needs.
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(_status)
