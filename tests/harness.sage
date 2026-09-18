r"""
Shared helpers for the test suite.

Loaded by ``tests/run_tests.sage`` and ``tests/make_references.sage``, both of
which are run **from the repository root**::

    sage tests/run_tests.sage

The scripts under test resolve their own paths relative to the current
directory, and they are not modules: they are ``load``ed, and a ``load``ed file
runs in the globals of whoever loads it. Three consequences shape this file.

- The example scripts, and with them ``code/sieve_integral.py``, are loaded
  into the globals of the test run itself. They have to be: the library sends
  the pieces of a balanced polytope to a ``multiprocessing.Pool``, which
  pickles its worker by name, and the name is only found if it lives in
  ``__main__``. Loading the library into a namespace of its own makes every
  example of more than twenty pieces fail with ``PicklingError``.
- ``Computation of h/computation_of_h.sage`` is loaded into a namespace of its
  own instead, because it would otherwise overwrite the library: both define
  ``are_inequalities_compatible``, with different arguments. Nothing the three
  h anchors reach uses the pool, so the isolation costs nothing there.
- Values are computed through ``example_value``, which caches them, so that
  several tests can look at the same example without paying for it twice.

Because the example scripts share the globals of the test run, nothing here or
in a ``test_*.sage`` file may bind a name they use -- in particular the
symbolic variables ``t``, ``u``, ``v``, ``w``, ``z``, ``x1`` ... ``a7``, the
dict ``a`` of the Stadlmann script, or their ``test_*`` functions.
"""
import inspect
import json
import os
import shutil
import subprocess
from contextlib import contextmanager

from sage.repl.load import load as sage_load

EXAMPLES_DIR = "Examples"
H_DIR = "Computation of h"
CODE_DIR = "code"
REFERENCES = "tests/reference_values.json"

if not all(os.path.isdir(d) for d in (EXAMPLES_DIR, H_DIR, CODE_DIR)):
    raise RuntimeError("run the test suite from the repository root: "
                       f"{os.getcwd()} has no {EXAMPLES_DIR}/ or {H_DIR}/")


@contextmanager
def in_directory(path):
    """Run the body with ``path`` as the current directory."""
    here = os.getcwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(here)


def load_isolated(directory, filename, **seed):
    """
    Load ``filename`` from ``directory`` into a namespace of its own, and
    return that namespace. ``seed`` is bound before loading, which is how
    ``EXAMPLES_AUTORUN`` is set to False.
    """
    namespace = {}
    exec("from sage.all import *", namespace)
    namespace.update(seed)
    with in_directory(directory):
        sage_load(filename, namespace)
    return namespace


EXAMPLE_FILES = ["tests_chen.sage", "tests_ford_maynard.sage",
                 "tests_maynard.sage", "tests_merikoski.sage",
                 "tests_stadlmann.sage"]

_namespaces = {}


def examples_namespace():
    """
    The globals of this test run, with the example scripts and the library
    loaded into them. ``run_examples.sage`` itself is deliberately not loaded:
    it would bring its own ``run``, ``select`` and ``main`` into these globals,
    on top of the runner's.
    """
    if "examples" not in _namespaces:
        namespace = globals()
        namespace["EXAMPLES_AUTORUN"] = False
        with in_directory(EXAMPLES_DIR):
            for filename in EXAMPLE_FILES:
                load(filename)
        namespace["ALL_EXAMPLES"] = {}
        for registry in ("CHEN_EXAMPLES", "FORD_MAYNARD_EXAMPLES",
                         "MAYNARD_EXAMPLES", "MERIKOSKI_EXAMPLES",
                         "STADLMANN_EXAMPLES"):
            namespace["ALL_EXAMPLES"].update(namespace[registry])
        _namespaces["examples"] = namespace
    return _namespaces["examples"]


def h_namespace():
    """The namespace of ``Computation of h/computation_of_h.sage``."""
    if "h" not in _namespaces:
        _namespaces["h"] = load_isolated(
            H_DIR, "computation_of_h.sage", EXAMPLES_AUTORUN=False)
    return _namespaces["h"]


def library_namespace():
    """
    Where ``sieve_integral`` and friends live: the same globals as the example
    scripts, which load the library themselves.
    """
    return examples_namespace()


def example_names():
    """Every example name, the sieve ones first, then the three h anchors."""
    return list(examples_namespace()["ALL_EXAMPLES"]) + list(h_namespace()["H_EXAMPLES"])


def is_h_example(name):
    return name.startswith("h/")


_values = {}


def clear_values():
    """
    Forget the computed values, so that the next test does its own work.

    The cache is what makes the suite quick -- a regression test and a bound
    test on the same example pay for it once -- but it would make a timing run
    measure the cache rather than the computation, so the runner clears it
    before every test it times.
    """
    _values.clear()


def example_value(name, extra_bits=0):
    """
    The value of the named example, computed at most once per precision.

    ``extra_bits`` is added to the example's own default precision -- to both
    components of their ``(p0, p1)`` pair for the h anchors, which is why it is
    a number of extra bits rather than a precision: the two chains do not count
    precision the same way. The h examples return the five components of
    ``h_bar``; the value kept here is the fourth, the one the notebook reports.
    """
    key = (name, extra_bits)
    if key in _values:
        return _values[key]
    default = default_precision(name)
    if is_h_example(name):
        test = h_namespace()["H_EXAMPLES"][name]
        precision = tuple(p + extra_bits for p in default)
        # two of the three anchors print their progress by default
        with in_directory(H_DIR):
            value = test(precision=precision, verbose=False)[3]
    else:
        test = examples_namespace()["ALL_EXAMPLES"][name]
        with in_directory(EXAMPLES_DIR):
            value = test(precision=default + extra_bits)
    _values[key] = value
    return value


def default_precision(name):
    """The precision the example runs at when it is not given one."""
    if is_h_example(name):
        func = h_namespace()["H_EXAMPLES"][name]
    else:
        test = examples_namespace()["ALL_EXAMPLES"][name]
        # the Stadlmann examples are partial(test_stadlmann, j)
        func = getattr(test, "func", test)
    return inspect.signature(func).parameters["precision"].default


class SerialPool:
    """
    A stand-in for ``multiprocessing.Pool`` that maps in this process, as
    ``Examples/run_examples.sage`` has; repeated here so that the h namespace,
    which does not load that file, can be made serial too.
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
    """Run the body with every loaded namespace's ``pool`` mapping in process."""
    saved = []
    for namespace in _namespaces.values():
        if "pool" in namespace:
            saved.append((namespace, namespace["pool"]))
            namespace["pool"] = SerialPool
    try:
        yield
    finally:
        for namespace, pool in saved:
            namespace["pool"] = pool


def ball_to_json(x):
    """A real ball as exact rational strings, as in hbar-echantillon.json."""
    return {"value": str(x),
            "mid": str(x.mid().exact_rational()),
            "rad": str(x.rad().exact_rational()),
            "prec": int(x.parent().precision())}


def ball_from_json(entry):
    """
    The ball stored by ``ball_to_json``. ``add_error`` rounds the radius
    outward, so this contains what was stored rather than repeating it.
    """
    return RealBallField(entry["prec"])(QQ(entry["mid"])).add_error(QQ(entry["rad"]))


def read_references(path=REFERENCES):
    with open(path) as f:
        return json.load(f)


def provenance():
    """Commit, date and versions, recorded with every generated file."""
    def git(*command):
        return subprocess.run(("git",) + command, capture_output=True,
                              text=True).stdout.strip()
    from datetime import datetime
    import sage.version
    return {"commit": git("rev-parse", "--short", "HEAD"),
            "branch": git("rev-parse", "--abbrev-ref", "HEAD"),
            "date": datetime.now().astimezone().isoformat(timespec="seconds"),
            "sage": sage.version.version,
            "latte": shutil.which("integrate") or "not found",
            "cores": os.cpu_count()}
