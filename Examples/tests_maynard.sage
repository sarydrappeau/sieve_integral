r"""
Tests on the integrals I_1, ..., I_9 from J. Maynard, *Primes with restricted
digits* (doi:10.1007/s00222-019-00865-6).

Extracted from ``Tests-Maynard.ipynb``. Run from this directory, as the
library is loaded by relative path::

    cd Examples && sage tests_maynard.sage

or, from a Sage session started here::

    sage: load("tests_maynard.sage")
    sage: test_M_I1(precision=20)

The examples are also registered in ``MAYNARD_EXAMPLES``, which
``run_examples.sage`` collects. The default precisions are the ones the
notebook times each integral at.
"""

load_attach_path("../code")
load("sieve_integral.py")
from itertools import product as iterproduct

var("u, v, w, t, z")
th2 = 17/40
th1 = 9/25


def test_M_I1(precision=30, verbose=0):
    P = {th2 - th1 < v,
         v < u,
         u < th1,
         v < (1-u)/2,
         1 - th1 < u + v,
         u + v + w == 1}
    return sieve_integral_harman(P, (w, v), precision=precision,
                                 verbose=verbose)


def test_M_I2(precision=30, verbose=0):
    P = {th2 - th1 < v,
         v < u,
         u < th1,
         th2 < u + v,
         u + v < 1 - th2,
         1 - th1 < 2*v + u,
         2*v + u < 1,
         u + v + w == 1}
    return sieve_integral(P, precision=precision, verbose=verbose)


def test_M_I3(precision=25, verbose=0):
    P = {th2 - th1 < v,
         v < u,
         u < th1,
         th2 < u + v,
         u + v < 1 - th2,
         1 - th1 < u + 2*v,
         u + 2*v < 1,
         v < w,
         w < (1-u-v)/2,
         v + w < th1,
         u + v + w + t == 1}
    return sieve_integral_harman(P, (t, w), precision=precision,
                                 verbose=verbose)


def test_M_I4(precision=22, verbose=0):
    P = {th2 - th1 < v,
         v < u,
         u < th1,
         th2 < u + v,
         u + v < 1 - th2,
         1 - th1 < u + 2*v,
         u + 2*v < 1,
         v < w,
         w < (1-u-v)/2,
         v + w > th2,
         u + v + w + t == 1}
    return sieve_integral_harman(P, (t, w), precision=precision,
                                 verbose=verbose)


def test_M_I5(precision=20, verbose=0):
    P_base = [th2 - th1 < t,
              t < w,
              w < v,
              v < u,
              u < th1,
              u + 2*v < 1 - th1,
              u + v + 2*w < 1,
              u + v + w + 2*t < 1,
              th2 < u + v < 1 - th2,
              u + v + w + t + z == 1]
    Plist = []
    vs = [u, v, w, t]
    sums = [vs[i] + vs[j] for i in range(4) for j in range(i+1, 4)]
    for eps in iterproduct([0, 1], repeat=6):
        pile = []
        for i in range(6):
            if eps[i] == 0:
                pile.append(sums[i] < th1)
            else:
                pile.append(sums[i] > th2)
        Plist.append(P_base + pile)

    res = 0
    for Q in Plist:
        res += sieve_integral_harman(Q, (z, t), precision=precision,
                                     verbose=verbose)
    return res


def test_M_I6(precision=20, verbose=0):
    I6_base = [th2 - th1 < t,
               t < w,
               w < v,
               v < u,
               u < th1,
               u + v < th1,
               u + v + w + t + z == 1]
    Plist = []
    ieqs_sum4 = [[u + v + w + t < th1],
                 [th2 < u + v + w + t, u + v + w + t < 1 - th2],
                 [u + v + w + t > 1 - th1]]
    vs = [u, v, w, t]
    sums = [u + v + w + t - x for x in vs]
    for ieqs_sum in ieqs_sum4:
        for eps in iterproduct([0, 1], repeat=4):
            pile = []
            for i in range(4):
                if eps[i] == 0:
                    pile.append(sums[i] < th1)
                else:
                    pile.append(sums[i] > th2)
            Plist.append(I6_base + ieqs_sum + pile)

    res = 0
    for Q in Plist:
        res += sieve_integral_harman(Q, (z, t), precision=precision,
                                     verbose=verbose)
    return res


def test_M_I7(precision=20, verbose=0):
    P = {th2 < u,
         u < 1/2,
         th2 - th1 < v,
         v < (1-u)/2,
         1 - th1 < u + v,
         u + v + w == 1}
    return sieve_integral_harman(P, (w, v), precision=precision,
                                 verbose=verbose)


def test_M_I8(precision=20, verbose=0):
    P = {th2 < u,
         u < 1/2,
         th2 - th1 < v,
         v < (1-u)/2,
         th2 < u + v,
         u + v < 1 - th2,
         1 - th1 < u + 2*v,
         u + v + w == 1}
    return sieve_integral_harman(P, (w, v), precision=precision,
                                 verbose=verbose)


def test_M_I9(precision=20, verbose=0):
    P_base = [th2 - th1 < t,
              t < w,
              w < v,
              th2 < u,
              u < 1/2,
              u + 2*v < 1 - th1,
              u + v + 2*w < 1,
              u + v + w + 2*t < 1,
              th2 < u + v,
              u + v < 1 - th2,
              u + v + w + t + z == 1]
    Plist = []
    vs = [u, v, w, t]
    sums = [vs[i] + vs[j] for i in range(4) for j in range(i+1, 4)]
    for eps in iterproduct([0, 1], repeat=6):
        pile = []
        for i in range(6):
            if eps[i] == 0:
                pile.append(sums[i] < th1)
            else:
                pile.append(sums[i] > th2)
        Plist.append(P_base + pile)

    res = 0
    for Q in Plist:
        res += sieve_integral_harman(Q, (z, t), precision=precision,
                                     verbose=verbose)
    return res


MAYNARD_EXAMPLES = {
    "maynard/I1": test_M_I1,
    "maynard/I2": test_M_I2,
    "maynard/I3": test_M_I3,
    "maynard/I4": test_M_I4,
    "maynard/I5": test_M_I5,
    "maynard/I6": test_M_I6,
    "maynard/I7": test_M_I7,
    "maynard/I8": test_M_I8,
    "maynard/I9": test_M_I9,
}


# See the note at the end of tests_chen.sage on this guard. Note that this
# file binds the symbolic variables u, v, w, t and z in the loading namespace:
# anything that loads it must not rebind those names.
import os as _os
import sys as _sys

if (_os.path.basename(_sys.argv[0]).startswith("tests_maynard")
        and globals().get("EXAMPLES_AUTORUN", True)):
    for _name, _test in MAYNARD_EXAMPLES.items():
        print(_name, "=", _test())
