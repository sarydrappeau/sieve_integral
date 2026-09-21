r"""
Tests from Section 2.5.1 of the preprint -- Chen's work.

Extracted from ``Tests-Chen.ipynb``. Run from this directory, as the library
is loaded by relative path::

    cd Examples && sage tests_chen.sage

or, from a Sage session started here::

    sage: load("tests_chen.sage")
    sage: test_chen_C1(precision=40)

The examples are also registered in ``CHEN_EXAMPLES``, which
``run_examples.sage`` collects; see that file for the timing harness. The
default precisions are the ones the notebook uses.
"""

load_attach_path("../code")
load("sieve_integral.py")

var("t1, t2, t3")


# Integral C1 = I^*(P), where P = {1/10 < t1 < 1/3 < t2 < t3}.

chen_P0 = {t1 > 1/10,
           t1 < 1/3,
           t2 > 1/3,
           t3 > t2,
           t1 + t2 + t3 == 1}


def test_chen_C1(precision=40, verbose=0):
    return sieve_integral(chen_P0, precision=precision, verbose=verbose)


# Integral C2 = I(P1) - I(P2)/2 - I(P3)/2, where
#   P1 = {1 < t3 < t2 - 1 < t1 - 2 < 2},
#   P2 = {1 < t3 < t2 - 1 < 3 - 10 t1 < 2},
#   P3 = {1 < t3 < t2 - 1 < 10 t1 - 2 < 2}.

chen_P1 = {1 < t3,
           t3 < t2 - 1,
           t2 - 1 < t1 - 2,
           t1 - 2 < 2}
chen_P2 = {1 < t3,
           t3 < t2 - 1,
           t2 - 1 < 3 - 10*t1,
           3 - 10*t1 < 2}
chen_P3 = {1 < t3,
           t3 < t2 - 1,
           t2 - 1 < 10*t1 - 2,
           10*t1 - 2 < 2}


def test_chen_C2(precision=40, verbose=0):
    v = [sieve_integral(Q, precision=precision, verbose=verbose)
         for Q in [chen_P1, chen_P2, chen_P3]]
    return v[0] - 1/2 * (v[1] + v[2])


CHEN_EXAMPLES = {
    "chen/C1": test_chen_C1,
    "chen/C2": test_chen_C2,
}


# Running this file directly computes every example above. ``__name__`` is
# "__main__" both when a .sage file is run as a script and when it is
# ``load``ed into a session, so it cannot tell the two apart; ``sys.argv[0]``
# names the preparsed script only in the first case. ``run_examples.sage``
# sets ``EXAMPLES_AUTORUN`` to False before loading this file, so that
# collecting the registry costs nothing.
import os as _os
import sys as _sys

if (_os.path.basename(_sys.argv[0]).startswith("tests_chen")
        and globals().get("EXAMPLES_AUTORUN", True)):
    for _name, _test in CHEN_EXAMPLES.items():
        print(_name, "=", _test())
