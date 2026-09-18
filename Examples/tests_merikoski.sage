r"""
Tests on the integrals F_1, ..., F_6 from J. Merikoski, *On the largest prime
factor of n^2 + 1* (doi:10.4171/JEMS/1216).

Extracted from ``Tests-Merikoski.ipynb``. Run from this directory, as the
library is loaded by relative path::

    cd Examples && sage tests_merikoski.sage

or, from a Sage session started here::

    sage: load("tests_merikoski.sage")
    sage: test_JM_F1(precision=20)

Compared with Merikoski's notation we change variables

    (u, v) = (alpha - beta, beta)   <=>   (alpha, beta) = (u + v, v)

so that the argument of the Buchstab function becomes u/v, which is the format
``sieve_integral_harman`` expects. In F4 the same is done with

    beta_i = v_i,   alpha - beta_1 - beta_2 - beta_3 = u.

The notebook has no F5. The examples are also registered in
``MERIKOSKI_EXAMPLES``, which ``run_examples.sage`` collects; the default
precisions are the ones the notebook uses.
"""

load_attach_path("../code")
load("sieve_integral.py")
from itertools import product as iterproduct

var("u, v, v1, v2, v3")


def test_JM_F1(precision=50, verbose=0):
    alpha = u + v
    beta = v
    sigma = (2-alpha)/3
    xi = 3/2 - alpha

    P11 = {1 < alpha,
           alpha < 17/16,
           sigma < beta,
           beta < alpha - 2*sigma}
    P12 = {1 < alpha,
           alpha < 17/16,
           xi < beta,
           beta < alpha / 2}

    factor2 = alpha

    F11 = sieve_integral_harman(P11, (u, v), factor2, precision=precision,
                                verbose=verbose)
    F12 = sieve_integral_harman(P12, (u, v), factor2, precision=precision,
                                verbose=verbose)
    return F11 + F12


def test_JM_F2(precision=50, verbose=0):
    alpha = u + v
    beta = v
    sigma = (2-alpha)/3

    P2 = {17/16 < alpha,
          alpha < 8/7,
          sigma < beta,
          beta < alpha/2}

    factor2 = alpha

    return sieve_integral_harman(P2, (u, v), factor2, precision=precision,
                                 verbose=verbose)


def test_JM_F3(precision=50, verbose=0):
    alpha = u + v
    beta = v
    sigma = (2-alpha)/3

    P3 = {8/7 < alpha,
          alpha < 7/6,
          sigma < beta,
          beta < alpha/2}

    factor2 = alpha

    return sieve_integral_harman(P3, (u, v), factor2, precision=precision,
                                 verbose=verbose)


def test_JM_F4(precision=50, verbose=0):
    alpha = u + v1 + v2 + v3
    beta1 = v1
    beta2 = v2
    beta3 = v3

    sigma = (2-alpha)/3
    gamma = (2-alpha)/3 - alpha + 1

    P4_base = [8/7 < alpha,
               alpha < 7/6,
               gamma < beta3,
               beta3 < beta2,
               beta2 < beta1,
               beta1 < alpha - 1]

    Plist = []
    for eps in iterproduct([False, True], repeat=4):
        Plist.append(P4_base +
                     ([beta1 + beta2 < alpha - 1] if eps[0]
                      else [beta1 + beta2 > sigma]) +
                     ([beta1 + beta3 < alpha - 1] if eps[1]
                      else [beta1 + beta3 > sigma]) +
                     ([beta2 + beta3 < alpha - 1] if eps[2]
                      else [beta2 + beta3 > sigma]) +
                     ([beta1 + beta2 + beta3 < alpha - 1] if eps[3]
                      else [beta1 + beta2 + beta3 > sigma]))
    factor2 = alpha

    res = 0
    for Q in Plist:
        res += sieve_integral_harman(Q, (u, v3), factor2,
                                     precision=precision, verbose=verbose)
    return res


def test_JM_F6(precision=30, verbose=0):
    alpha = u + v
    beta = v
    sigma = (2-alpha)/3

    P6 = {7/6 < alpha,
          alpha < 5/4,
          alpha - 1 < beta,
          beta < sigma}

    factor2 = alpha

    return sieve_integral_harman(P6, (u, v), factor2, precision=precision,
                                 verbose=verbose)


MERIKOSKI_EXAMPLES = {
    "merikoski/F1": test_JM_F1,
    "merikoski/F2": test_JM_F2,
    "merikoski/F3": test_JM_F3,
    "merikoski/F4": test_JM_F4,
    "merikoski/F6": test_JM_F6,
}


# See the note at the end of tests_chen.sage on this guard. Note that this
# file binds the symbolic variables u, v, v1, v2 and v3 in the loading
# namespace: anything that loads it must not rebind those names.
import os as _os
import sys as _sys

if (_os.path.basename(_sys.argv[0]).startswith("tests_merikoski")
        and globals().get("EXAMPLES_AUTORUN", True)):
    for _name, _test in MERIKOSKI_EXAMPLES.items():
        print(_name, "=", _test())
