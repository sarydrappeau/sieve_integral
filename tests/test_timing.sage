r"""
Micro-benchmarks for the parts later steps rewrite, which the published
examples do not isolate. They are tagged ``timing`` and are only run by

    sage tests/run_tests.sage --only-timing [--long] [--time OUT.json]
                              [--compare BASE.json]

Each one does a fixed amount of work and checks nothing beyond finishing; the
runner measures it. The polytopes have small denominators, as the floating
point LP of the library needs.

What each group is for:

- ``balanced-integrate``: the Taylor expansion and the LattE call of
  ``balanced_polytope_integrate`` on one balanced piece, at three dimensions
  and two precisions -- the target of the closed-form Taylor coefficients
  (step 10) and of the Chebyshev expansion (step 13).
- ``border-integrate``: the same on a piece cut out by ``sum x_i = 1``, where
  the last coordinate is eliminated -- the target of step 11.
- ``dde``: building the Buchstab and Dickman approximants up to u = 10, the
  target of step 14.
"""
import itertools


def balanced_box(dimension, ratio=QQ(13)/10):
    """
    A balanced box: the j-th coordinate runs from 1 + j/10 to ratio times
    that, so that max/min is the same in every coordinate, as the pieces
    ``balance_polytope`` produces.
    """
    intervals = [(1 + QQ(j) / 10, (1 + QQ(j) / 10) * ratio)
                 for j in range(dimension)]
    return Polyhedron(vertices=list(itertools.product(*intervals)),
                      base_ring=QQ)


def border_box(dimension, low=QQ(1)/5, high=QQ(13)/50):
    """
    A piece of the hyperplane sum x_i = 1, with the first ``dimension - 1``
    coordinates in [low, high]; the last one is what is left.
    """
    eqn = [-1] + [1] * dimension
    ieqs = []
    for j in range(dimension - 1):
        row = [0] * (dimension + 1)
        row[0], row[j + 1] = -low, 1
        ieqs.append(list(row))
        row = [0] * (dimension + 1)
        row[0], row[j + 1] = high, -1
        ieqs.append(list(row))
    polytope = Polyhedron(eqns=[eqn], ieqs=ieqs, base_ring=QQ)
    assert polytope.dimension() == dimension - 1, (
        f"the border box in dimension {dimension} came out with dimension "
        f"{polytope.dimension()}")
    return polytope


def integrate_piece(polytope, precision, degree=10):
    library = library_namespace()
    summary = library["PolytopeSummary"](polytope)
    return library["balanced_polytope_integrate"](
        summary, truncation_degree=degree,
        scalar_field=RealBallField(precision))


def balanced_integrate_benchmark(dimension, precision):
    polytope = balanced_box(dimension)

    def benchmark():
        integrate_piece(polytope, precision)
    return benchmark


def border_integrate_benchmark(dimension, precision):
    polytope = border_box(dimension)

    def benchmark():
        integrate_piece(polytope, precision)
    return benchmark


def dde_benchmark(name, precision, argument=10, repeats=20):
    """
    Build the solution and extend it to u = ``argument``, ``repeats`` times
    over: one construction takes a few milliseconds even at 200 bits, which is
    too close to the noise for a before-and-after comparison to mean anything.
    """
    def benchmark():
        library = library_namespace()
        for _ in range(repeats):
            solution = library[name](precision=precision)
            solution(argument)
        assert solution.approximants_list(), f"{name} has no approximants"
    return benchmark


TESTS = {}

for _dimension in (3, 4, 5):
    for _precision in (53, 200):
        TESTS[f"timing/balanced-integrate/dim{_dimension}-prec{_precision}"] = (
            balanced_integrate_benchmark(_dimension, _precision), ("timing",))

for _precision in (53, 200):
    TESTS[f"timing/border-integrate/dim4-prec{_precision}"] = (
        border_integrate_benchmark(4, _precision), ("timing",))

for _name in ("BuchstabB", "DickmanRho"):
    for _precision in (50, 100, 200):
        TESTS[f"timing/dde/{_name}-prec{_precision}"] = (
            dde_benchmark(_name, _precision), ("timing",))
