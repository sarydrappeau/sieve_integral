r"""
Tests on several integrals from J. Stadlmann, *On the mean square gap
between primes* (arXiv:2212.10867).

Extracted from ``Tests-Stadlmann.ipynb``. Run from this directory, as
the library is loaded by relative path::

    cd Examples && sage tests_stadlmann.sage

or, from a Sage session started here::

    sage: load("tests_stadlmann.sage")
    sage: test_stadlmann(0, precision=20)

The examples are also registered in ``STADLMANN_EXAMPLES``, which
``run_examples.sage`` collects.

``STADLMANN_POLYTOPES[j]`` is the j-th polytope and
``STADLMANN_BOUNDS[j]`` the upper bound obtained in the paper for the
corresponding integral -- which is the easiest way to locate a given
polytope in the paper. The variable ``a[i]`` corresponds to
`\alpha_i` there.

Note that for each polytope a dummy variable `a_{k+1} = 1 - \sum_{i\le k} a_i`
is added, to match the shape of the integrals ``sieve_integral_harman``
computes.

Indices 40 and 41 are the six-dimensional integrals on page 67; they take
substantially longer than the rest, and the notebook runs them apart.
"""

load_attach_path("../code")
load("sieve_integral.py")
from functools import partial

var("a1, a2, a3, a4, a5, a6, a7")

a = {1: a1,
     2: a2,
     3: a3,
     4: a4,
     5: a5,
     6: a6,
     7: a7}


def polytope_from_bounds(L, var=0):
    """
    Here L is a lists [c1, c2, ..., c_k]
    where k is the number of variables in the to-be-computed integral
    (full-dimensional). Each c_i is a list [m, M], where m and M are
    themselves a list m = [m1, ... m_r],  M = [M_1, ..., M_s],
    corresponding to conditions  max(m) <= a_i <= min(M)
    """
    P = {add(a[j] for j in range(1, var+2)) == 1}
    for j in range(len(L)):
        for u in L[j][0]:
            P.add(a[j+1] > u)
        for u in L[j][1]:
            P.add(a[j+1] < u)
    return P


def _build_stadlmann_polytopes():
    """
    Return the pair ``(polytope_list, bounds_from_arxiv)``.

    Kept as a function rather than run at load time, so that importing
    this file has no side effect beyond binding the two module-level
    lists below.
    """
    polytope_list = []
    bounds_from_arxiv = []

    def add_polytope(L, var=0):
        polytope_list.append(polytope_from_bounds(L, var=var))

    # 0
    #
    # For instance, the following example corresponds to the polytope
    # defined by
    #     0.345 <= a1 <= 0.375
    #     0.485 - a1 <= a2 <= (0.655 - a1)/2
    #     0.08 <= a3
    #     0.08 <= a4 <= min( a3, (1-a1-a2-a3)/2 )
    # This is the first integral on page 58.
    add_polytope([[[0.345], [0.375]],
                  [[0.485-a1], [(0.655-a1)/2]],
                  [[0.08], [a2]],
                  [[0.08], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.012)

    add_polytope([[[0.285], [0.315]],
                  [[0.08], [0.405-a1]],
                  [[0.08], [a2]],
                  [[0.08], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.003)

    add_polytope([[[0.655/3], [0.285]],
                  [[(0.655-a1)/2], [a1, 0.526-a1]]], var=2)
    bounds_from_arxiv.append(0.01)

    add_polytope([[[0.08], [0.285]],
                  [[0.08], [a1, (0.655-a1)/2]],
                  [[0.08], [a2]],
                  [[0.08, (0.655-a1-a2-a3)/2], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.296)

    add_polytope([[[0.455], [0.5]],
                  [[0.62-a1], [(1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.2029)

    add_polytope([[[0.455], [0.475]],
                  [[(0.685-a1)/2], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.0099)

    add_polytope([[[0.455], [0.5]],
                  [[0.075], [0.58-a1, (0.685-a1)/2]],
                  [[0.075], [a2]],
                  [[0.075], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0038)

    add_polytope([[[0.42], [0.455]],
                  [[0.685-a1], [(1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.0345)

    add_polytope([[[0.42], [0.455]],
                  [[(0.685-a1)/2], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.0502)

    add_polytope([[[0.42], [0.455]],
                  [[0.105], [(0.685-a1)/2]],
                  [[0.105], [a2]],
                  [[0.105], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0004)

    # 10
    add_polytope([[[0.315], [0.38]],
                  [[0.62-a1], [a1, (1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.0889)

    add_polytope([[[0.315], [0.38]],
                  [[0.145, (0.62-a1)/2], [0.545-a1]]], var=2)
    bounds_from_arxiv.append(0.1993)

    add_polytope([[[0.31], [0.315]], [[0.62-a1], [a1]]], var=2)
    bounds_from_arxiv.append(0.0007)

    add_polytope([[[0.29], [0.315]], [[(0.62-a1)/2], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.1343)

    add_polytope([[[0.29], [0.315]], [[0.42-a1], [(0.62-a1)/2]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0020)

    add_polytope([[[0.29], [0.305]], [[0.075], [(0.62-a1)/2]],
                  [[0.075], [a2]],
                  [[0.075], [0.38-a1, a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0092)

    add_polytope([[[0.075], [0.29]], [[0.455-a1], [a1, (0.685-a1)/2]],
                  [[0.075], [0.38-a1]],
                  [[0.075, (0.62-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0116)

    add_polytope([[[0.075], [0.29]], [[0.455-a1], [a1, (0.685-a1)/2]],
                  [[0.42-a1], [a2]],
                  [[0.075, (0.62-a1-a2-a3)/2],
                   [a3, 0.38-a1, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0129)

    add_polytope([[[0.075], [0.029]], [[0.455-a1], [a1, (0.685-a1)/2]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1, (0.62-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0027)

    add_polytope([[[0.105], [0.29]], [[0.42-a1], [a1, 0.455-a1]],
                  [[0.105], [a2]],
                  [[0.42-a1, (0.62-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0012)

    # 20
    add_polytope([[[0.075], [0.29]], [[0.075, 0.315-a1], [a1, 0.38-a1]],
                  [[0.075], [a2]],
                  [[0.455-a2-a3], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0048)

    add_polytope([[[0.075], [0.29]], [[0.075, 0.315-a1], [a1, 0.38-a1]],
                  [[0.075], [a2]],
                  [[0.455-a1-a2], [a3, 0.38-a2-a3]]], var=4)
    bounds_from_arxiv.append(0.0252)

    add_polytope([[[0.075], [0.29]], [[0.075], [a1, 0.315-a1]],
                  [[0.075], [a2]],
                  [[0.075, (0.62-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0127)

    add_polytope([[[0.435], [0.5]], [[0.635-a1], [(1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.2182)

    add_polytope([[[0.435], [0.5]],
                  [[0.07], [0.58-a1, 0.105, (0.67-a1)/2]],
                  [[0.07], [a2]],
                  [[0.07], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0083)

    add_polytope([[[0.33], [0.365]], [[0.635-a1], [a1, (1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.0367)

    add_polytope([[[0.33], [0.365]], [[0.1524], [0.565-a1]]], var=2)
    bounds_from_arxiv.append(0.1186)

    add_polytope([[[0.33], [0.365]], [[0.435-a1], [0.1524]],
                  [[0.435-a1], [a2]],
                  [[0.435-a1], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0211)

    add_polytope([[[0.305], [0.33]], [[0.635-a1], [a1]]], var=2)
    bounds_from_arxiv.append(0.0043)

    add_polytope([[[0.305], [0.33]], [[(0.635-a1)/2], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.1178)

    # 30
    add_polytope([[[0.305], [0.33]], [[0.42-a1], [(0.635-a1)/2]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0062)

    add_polytope([[[0.2099], [0.305]], [[0.2099], [a1, 0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.1723)

    add_polytope([[[0.21], [0.305]], [[0.42-a1], [0.2099]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0104)

    add_polytope([[[0.21], [0.305]], [[0.42-a1], [0.2099]],
                  [[0.42-a1], [a2]],
                  [[0.07], [a3, (1-a1-a2-a3)/2, 0.365-a1]]], var=4)
    bounds_from_arxiv.append(0.0212)

    add_polytope([[[0.21], [0.305]], [[0.42-a1], [0.2099]],
                  [[0.07], [a2, 0.365-a1]],
                  [[0.07], [a3, 0.365-a1]]], var=4)
    bounds_from_arxiv.append(0.0397)

    add_polytope([[[0.07], [0.295]], [[0.07], [a1, 0.365-a1]],
                  [[0.07], [a2]],
                  [[(0.635-a1-a2-a3)/2, 0.42-a2-a3],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0105)

    add_polytope([[[0.36], [0.5]], [[0.07], [a1, (0.71-a1)/2]],
                  [[0.07], [a2, 0.64-a1-a2]],
                  [[0.07], [a3]]], var=4)
    bounds_from_arxiv.append(0.079)

    add_polytope([[[0.07], [0.29]], [[(0.71-a1)/2], [a1]]], var=2)
    bounds_from_arxiv.append(0.08)

    add_polytope([[[0.22], [0.29]], [[0.07], [a1, (0.71-a1)/2]],
                  [[0.07], [a2]],
                  [[0.36-a1, (0.71-a1-a2-a3)/2], [a3]]], var=4)
    bounds_from_arxiv.append(0.112)

    add_polytope([[[0.07], [0.22]], [[0.07], [a1, (0.71-a1)/2]],
                  [[0.07], [a2]],
                  [[0.07, (0.71-a1-a2-a3)/2], [a3]]], var=4)
    bounds_from_arxiv.append(0.063)

    # 40
    add_polytope([[[0.18], [0.29]], [[0.145], [0.29, a1, (0.71-a1)/2]],
                  [[0.07], [0.145, a2]],
                  [[0.07], [a3, (0.71-a1-a2-a3)/2]],
                  [[0.07], [a4]],
                  [[0.07], [a5]]], var=6)
    bounds_from_arxiv.append(0.056)

    polytope_list.append({a6 > 0.07, a5 > a6, a4 > a5, a3 > a4, a2 > a3,
                          a1 > a2, 1/2 > a1,
                          a1 < 0.29,
                          a1 + 2*a2 < 0.71,
                          a1 + a2 + a3 + 2*a4 < 0.71,
                          a2 > 0.145,
                          a1 + a2 + a3 + a4 + a5 + a6 + a7 == 1})
    bounds_from_arxiv.append(0.035)

    add_polytope([[[0.475], [0.5]], [[0.6-a1], [(1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.166)

    add_polytope([[[0.3], [0.4]], [[0.6-a1], [a1, (1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.187)

    add_polytope([[[0.475/2], [0.385]],
                  [[0.14, 0.475-a1], [a1, 0.525-a1]]], var=2)
    bounds_from_arxiv.append(0.302)

    add_polytope([[[0.335], [0.4]], [[0.475-a1], [0.14, 0.525-a1]],
                  [[0.075], [a2]],
                  [[0.075], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.032)

    add_polytope([[[0.075], [0.3255]], [[0.075], [a1, 0.4-a1]],
                  [[0.075], [a2]],
                  [[(0.615-a1-a2-a3)/2, 0.475-a1-a2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.07)

    add_polytope([[[0.075], [0.325]], [[0.075], [a1, 0.4-a1]],
                  [[0.075], [a2]],
                  [[0.075, (0.615-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2, 0.4-a1-a2]]], var=4)
    bounds_from_arxiv.append(0.01)

    add_polytope([[[0.42], [0.5]], [[0.645-a1], [(1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.2194)

    add_polytope([[[0.42], [0.48]], [[0.1], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.1769)

    # 50
    add_polytope([[[0.42], [0.5]], [[0.065], [0.1, 0.58-a1]],
                  [[0.065], [a2]],
                  [[0.065], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0170)

    add_polytope([[[0.3225], [0.355]], [[0.645-a1], [a1, (1-a1)/2]]], var=2)
    bounds_from_arxiv.append(0.0191)

    add_polytope([[[0.325], [0.355]], [[(0.645-a1)/2], [0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.1266)

    add_polytope([[[0.325], [0.355]], [[0.42-a1], [(0.645-a1)/2]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1], [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0282)

    add_polytope([[[0.2099], [0.325]], [[0.2099], [a1, 0.58-a1]]], var=2)
    bounds_from_arxiv.append(0.2102)

    add_polytope([[[0.21], [0.325]], [[0.42-a1], [0.21]],
                  [[0.065], [a2, 0.355-a1]],
                  [[0.065, (0.42-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0249)

    add_polytope([[[0.21], [0.325]], [[0.42-a1], [0.21]],
                  [[0.42-a1], [a2]],
                  [[0.065, (0.42-a2-a3)/2],
                   [a3, 0.355-a1, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0191)

    add_polytope([[[0.21], [0.325]], [[0.42-.0], [0.21]],
                  [[0.42-a1], [a2]],
                  [[0.42-a1, (0.42-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0280)

    add_polytope([[[0.065], [0.325]], [[0.065], [a1, 0.355-a1]],
                  [[0.065], [a2]],
                  [[0.065, (0.645-a1-a2-a3)/2],
                   [a3, (1-a1-a2-a3)/2]]], var=4)
    bounds_from_arxiv.append(0.0180)

    # Three of the polytopes need further corrections (mismatching
    # between the p_i-sum conditions and the integral condition).

    # This one is the one whose integral is bounded by 0.0027
    polytope_list[18] = {a4 > 0.075,
                         a3 > a4,
                         a2 > a3,
                         a1 > a2,
                         1/2 > a1,
                         0.29 > a1,
                         0.2275 > a2,
                         a1 + 2*a2 < 0.685,
                         a1 + a2 + a3 + 2*a4 > 0.62,
                         a1 + a2 > 0.455,
                         a1 + a3 > 0.42,
                         a1 + a4 > 0.42,
                         a1 + a2 + a3 + a4 + a5 == 1}

    # This one is the one whose integral is bounded by 0.0180
    polytope_list[58] = {a4 > 0.065,
                         a3 > a4,
                         a2 > a3,
                         a1 > a2,
                         1/2 > a1,
                         a1 < 0.325,
                         a1 + a2 < 0.355,
                         a1 + a2 + a3 + 2*a4 > 0.645,
                         a1 + a2 + a3 + a4 + a5 == 1}

    # This one is the one whose integral is bounded by 0.0280
    polytope_list[57] = {a4 > 0.065,
                         a3 > a4,
                         a2 > a3,
                         a1 > a2,
                         1/2 > a1,
                         a1 < 0.325,
                         a2 < 0.2099,
                         a1 + a2 > 0.42,
                         a2 + a3 + 2*a4 < 0.42,
                         a1 + a3 > 0.42,
                         a1 + a4 > 0.42,
                         a1 + a2 + a3 + a4 + a5 == 1}

    return polytope_list, bounds_from_arxiv


STADLMANN_POLYTOPES, STADLMANN_BOUNDS = _build_stadlmann_polytopes()

# The six-dimensional integrals on page 67. They take substantially
# longer than the others, and the notebook runs them apart from the rest.
STADLMANN_SLOW = (40, 41)


def test_stadlmann(j, precision=20, verbose=0):
    """
    Compute the j-th Stadlmann integral, as ``work`` does in the
    notebook: the Buchstab variables are the last two coordinates.
    """
    Q = STADLMANN_POLYTOPES[j]
    eqns, ieqs, v = symbolic_to_eqns(Q)
    P = Polyhedron(eqns=eqns, ieqs=ieqs, base_ring=QQ)
    k = P.ambient_dimension()
    return sieve_integral_harman(Q, (a[k], a[k-1]), precision=precision,
                                 verbose=verbose)


STADLMANN_EXAMPLES = {
    f"stadlmann/{j:02d}": partial(test_stadlmann, j)
    for j in range(len(STADLMANN_POLYTOPES))
}


# See the note at the end of tests_chen.sage on this guard. Running this file
# directly reproduces the notebook's own loop: every integral but the two slow
# ones, against the bound the paper gives for it. The names bound below are
# prefixed with an underscore because the test files put their symbolic
# variables in the loading namespace, and t1, t2 are Chen's.
import os as _os
import sys as _sys

if (_os.path.basename(_sys.argv[0]).startswith("tests_stadlmann")
        and globals().get("EXAMPLES_AUTORUN", True)):
    from time import time as _time

    _K = RealBallField(precision=20)
    _computation_time = {}
    for _j in range(len(STADLMANN_POLYTOPES)):
        if _j in STADLMANN_SLOW:
            # Run apart; they take quite longer.
            continue
        _start = _time()
        _value = _K(test_stadlmann(_j))
        _computation_time[_j] = round(_time() - _start, 2)
        print(f" I[{_j}] <= {_value} versus {STADLMANN_BOUNDS[_j]} "
              f"(time = {_computation_time[_j]})")
        print(f"  difference {STADLMANN_BOUNDS[_j] - _value.upper()}\n")
    print("total time", add(_computation_time.values()))
