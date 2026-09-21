r"""
Tests on five integrals from K. Ford and J. Maynard, *On the theory of prime
producing sieves* (arXiv:2407.14368).

Extracted from ``Tests-Ford-Maynard.ipynb``. Run from this directory, as the
library is loaded by relative path::

    cd Examples && sage tests_ford_maynard.sage

or, from a Sage session started here::

    sage: load("tests_ford_maynard.sage")
    sage: test_FM_I3(precision=20)

The examples are also registered in ``FORD_MAYNARD_EXAMPLES``, which
``run_examples.sage`` collects. The default precisions are the ones the
notebook times each integral at.
"""

load_attach_path("../code")
load("sieve_integral.py")

# As in the notebook: a RealNumber run through QQ, that is the continued
# fraction reconstruction of 0.16623, not its exact binary value.
nu = QQ(0.16623)

var("x1, x2, x3, x4, x5, x6")


def test_FM_I3(precision=33, verbose=0):
    P = {nu < x1,
         x1 < x2,
         x2 < x3,
         x3 < 1/2,
         x1 + x2 + x3 == 1}
    I3a = sieve_integral(P, precision=precision, verbose=verbose,
                         facteur=1.15)

    P = {nu < x1,
         x1 < x2,
         x2 < x3,
         x1 + x2 < (1-nu)/2,
         x1 + x2 + x3 == 1}
    I3b = sieve_integral(P, precision=precision, verbose=verbose,
                         facteur=1.15)

    return -2*I3a - I3b


FM_I4_P = [{nu < x1,
            x1 < x2,
            x2 < x3,
            x3 < x4,
            x1 + x2 + x3 + x4 == 1,
            s < (1-nu)/2} for s in [x1+x2, x1+x3, x2+x3]]


def test_FM_I4(precision=34, verbose=0):
    return -add([sieve_integral(Q, precision=precision, verbose=verbose)
                 for Q in FM_I4_P])


FM_I5_base = [nu < x1,
              x1 < x2,
              x2 < x3,
              x3 < x4,
              x4 < x5,
              x1 + x2 + x3 + x4 + x5 == 1]

FM_I5_P = ([FM_I5_base]
           + [FM_I5_base + [s < (1-nu)/2]
              for s in (x3+x4, x1+x5, x2+x5, x3+x5, x4+x5)])


def test_FM_I5(precision=31, verbose=0):
    res = [sieve_integral(Q, precision=precision, verbose=verbose)
           for Q in FM_I5_P]
    return res[0] - add(res[1:])


def test_FM_I5b(precision=37, verbose=0):
    """This one is for I_5'."""
    tx = [x1, x2, x3, x4, x5]
    I5b_P = [FM_I5_base + [tx[i] + tx[j] > 1/2]
             for i in range(5) for j in range(i+1, 5)]
    I5b = 0
    for Q in I5b_P:
        I5b += -2 * sieve_integral(Q, precision=precision, verbose=verbose)
    return I5b


FM_tx6 = [x1, x2, x3, x4, x5, x6]

FM_I6_P = ([nu < x1]
           + [FM_tx6[i] < FM_tx6[i+1] for i in range(5)]
           + [sum(FM_tx6) == 1])


def test_FM_I6(precision=60, verbose=0):
    return -15 * sieve_integral(FM_I6_P, precision=precision, verbose=verbose)


FORD_MAYNARD_EXAMPLES = {
    "ford-maynard/I3": test_FM_I3,
    "ford-maynard/I4": test_FM_I4,
    "ford-maynard/I5": test_FM_I5,
    "ford-maynard/I5b": test_FM_I5b,
    "ford-maynard/I6": test_FM_I6,
}


# See the note at the end of tests_chen.sage on this guard.
import os as _os
import sys as _sys

if (_os.path.basename(_sys.argv[0]).startswith("tests_ford_maynard")
        and globals().get("EXAMPLES_AUTORUN", True)):
    for _name, _test in FORD_MAYNARD_EXAMPLES.items():
        print(_name, "=", _test())
