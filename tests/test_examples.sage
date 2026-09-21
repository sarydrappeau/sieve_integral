r"""
Tests on the published examples: ``Examples/run_examples.sage`` (80 of them)
and the three h anchors of ``Computation of h/computation_of_h.sage``.

Three kinds of check, in increasing order of what they would catch:

1. **regression** -- the value overlaps the one stored in
   ``tests/reference_values.json``. Those references are what the code gave on
   the branch that generated them; they are not certified values, and
   ``tests/make_references.sage`` regenerates them, which is an explicit act
   whose commit has to justify every entry that moved.
2. **self-consistency** -- the value overlaps the value at ten more bits of
   precision. Two enclosures of the same integral must meet, so this catches an
   error term that is accounted for too optimistically, which a regression test
   cannot.
3. **independent values** -- the closed forms of the library's docstrings, the
   outputs the h notebook has stored, and the upper bounds Stadlmann's paper
   gives for its polytopes.

The acceptance criterion throughout is that the result is an enclosure, never
that it repeats a previous value bit for bit.
"""
from fractions import Fraction

# Fast examples, one or more per source notebook, covering both the plain and
# the Buchstab chains; the rest are tagged long. Times from
# "Work plan/runs/02-examples-as-scripts/results-87f4475.json".
QUICK_EXAMPLES = [
    "chen/C1",
    "ford-maynard/I3", "ford-maynard/I5b", "ford-maynard/I6",
    "maynard/I1", "maynard/I2", "maynard/I3", "maynard/I7", "maynard/I8",
    "merikoski/F1", "merikoski/F6",
    "stadlmann/02", "stadlmann/04", "stadlmann/05", "stadlmann/09",
    "stadlmann/12", "stadlmann/25", "stadlmann/31", "stadlmann/37",
    "stadlmann/42", "stadlmann/51", "stadlmann/57",
    "h/0.4-0.6", "h/1over6-1over2",
]

# h/0.45-0.5 takes about eight seconds at its own precision and much longer at
# ten bits more, and the two six-dimensional Stadlmann integrals are slow
# enough that the notebook runs them apart: they keep their regression test but
# not the self-consistency one.
NO_SELF_CONSISTENCY = ["h/0.45-0.5", "stadlmann/40", "stadlmann/41"]

REFERENCE_DATA = read_references()


def regression_test(name):
    def test():
        entry = REFERENCE_DATA["examples"].get(name)
        assert entry is not None, (
            f"{name} has no reference value; run tests/make_references.sage")
        stored = ball_from_json(entry)
        value = example_value(name)
        assert value.overlaps(stored), (
            f"{name} = {value} does not overlap the reference {stored} "
            f"(from commit {REFERENCE_DATA['provenance']['commit']})")
    return test


def self_consistency_test(name):
    def test():
        coarse = example_value(name)
        fine = example_value(name, extra_bits=10)
        assert coarse.overlaps(fine), (
            f"{name} = {coarse} at its own precision does not overlap "
            f"{fine} at ten bits more: one of the two error terms is "
            f"accounted for too optimistically")
    return test


def stadlmann_bound_test(j):
    def test():
        name = f"stadlmann/{j:02d}"
        bound = examples_namespace()["STADLMANN_BOUNDS"][j]
        value = example_value(name)
        assert value.upper() <= bound, (
            f"{name} = {value} is not below the bound {bound} of the paper")
    return test


# ---------------------------------------------------------------- closed forms
#
# The values the docstrings of code/sieve_integral.py give, recomputed here as
# balls rather than copied as decimals: Arb integrates rigorously, so these are
# genuine enclosures of the true values.

def closed_form_log2():
    """Polyhedron([(1,), (2,)]): the integral is log 2."""
    library = library_namespace()
    with in_directory(CODE_DIR):
        value = library["sieve_integral"](Polyhedron([(1,), (2,)]), precision=50)
    true = RealBallField(200)(2).log()
    assert value.overlaps(true), f"{value} does not contain log 2 = {true}"


def closed_form_cube():
    """The unit cube at (2, 2, 2): the integral is log(3/2)^3."""
    library = library_namespace()
    polytope = vector((2, 2, 2)) + polytopes.hypercube(3, intervals='zero_one')
    with in_directory(CODE_DIR):
        value = library["sieve_integral"](polytope, precision=50)
    true = RealBallField(200)(3/2).log()**3
    assert value.overlaps(true), f"{value} does not contain log(3/2)^3 = {true}"


def closed_form_border():
    """
    The border polytope of the docstring: the integral is
    (1/2) * int_1^5 log(6 - x)/x dx.
    """
    library = library_namespace()
    polytope = Polyhedron(eqns=[(-1, 1, 1, 1)],
                          ieqs=[(-1/7, 1, 0, 0), (0, -1, 1, 0), (0, 0, -1, 1)])
    with in_directory(CODE_DIR):
        value = library["sieve_integral"](polytope, precision=50)
    true = (ComplexBallField(200).integral(lambda x, _: (6 - x).log() / x, 1, 5) / 2).real()
    assert value.overlaps(true), f"{value} does not contain {true}"


def closed_form_harman():
    """
    The docstring example of sieve_integral_harman, whose value is
    int_{1/3}^{1/2} (1+x)/(x(1-x)) dx + int_{1/4}^{1/3} (1+log(1/x-2))(1+x)/(x(1-x)) dx.
    """
    library = library_namespace()
    x, y = var("x, y")
    with in_directory(CODE_DIR):
        value = library["sieve_integral_harman"]({x + y == 1, 1/4 < x, x < y},
                                                 (y, x), 1 + x, precision=50)
    C = ComplexBallField(200)
    true = (C.integral(lambda t, _: (1 + t) / (t * (1 - t)), 1/3, 1/2)
            + C.integral(lambda t, _: (1 + (1/t - 2).log()) * (1 + t) / (t * (1 - t)),
                         1/4, 1/3)).real()
    assert value.overlaps(true), f"{value} does not contain {true}"


def closed_form_grid_rounding():
    """
    A two-dimensional box whose integral is log(285605/100000) * log(6/5).

    The ratio 285605/100000 lies just above the fourth power of 1.3 rounded to
    three decimals, and below the exact power, so the grid of
    ``balance_polytope`` has to reach it: otherwise a sliver of the polytope is
    dropped and the returned ball misses the true value.
    """
    library = library_namespace()
    polytope = Polyhedron(ieqs=[(-1, 1, 0), (QQ(285605/100000), -1, 0),
                                (-1, 0, 1), (QQ(6/5), 0, -1)])
    with in_directory(CODE_DIR):
        value = library["sieve_integral"](polytope, precision=30)
    R = RealBallField(200)
    true = R(285605/100000).log() * R(6/5).log()
    assert value.overlaps(true), f"{value} does not contain {true}"


# --------------------------------------------------------- the notebook's h
#
# The values printed in the outputs of "Computation of h(alpha, beta).ipynb",
# as the balls they print: our value has to meet them.

NOTEBOOK_H = {
    "h/0.4-0.6": ("0.511", "3.19e-4", 14),
    "h/1over6-1over2": ("0.182322", "6.90e-7", 20),
    "h/0.45-0.5": ("0.6", "0.0355", 14),
}


def notebook_h_test(name):
    def test():
        # QQ does not read a decimal string; Fraction does, exponent included.
        mid, rad, prec = NOTEBOOK_H[name]
        stored = RealBallField(prec)(QQ(Fraction(mid))).add_error(QQ(Fraction(rad)))
        value = example_value(name)
        assert value.overlaps(stored), (
            f"{name} = {value} does not overlap the notebook's {stored}")
    return test


# ------------------------------------------------------------------ registry

TESTS = {}

for _name in example_names():
    _tags = ("quick",) if _name in QUICK_EXAMPLES else ("long",)
    TESTS[f"regression/{_name}"] = (regression_test(_name), _tags)
    if _name not in NO_SELF_CONSISTENCY:
        TESTS[f"selfconsistency/{_name}"] = (self_consistency_test(_name), ("long",))

for _j in range(len(examples_namespace()["STADLMANN_BOUNDS"])):
    _tags = ("quick",) if f"stadlmann/{_j:02d}" in QUICK_EXAMPLES else ("long",)
    TESTS[f"bound/stadlmann/{_j:02d}"] = (stadlmann_bound_test(_j), _tags)

for _name in NOTEBOOK_H:
    _tags = ("quick",) if _name in QUICK_EXAMPLES else ("long",)
    TESTS[f"notebook/{_name}"] = (notebook_h_test(_name), _tags)

TESTS.update({
    "closedform/log2": (closed_form_log2, ("quick",)),
    "closedform/cube": (closed_form_cube, ("quick",)),
    "closedform/border": (closed_form_border, ("quick",)),
    "closedform/harman": (closed_form_harman, ("quick",)),
    "closedform/grid-rounding": (closed_form_grid_rounding, ("quick",)),
})


KNOWN_ISSUES = {}
