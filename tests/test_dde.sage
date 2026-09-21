r"""
Tests on ``code/number_theoretic_dde_solutions.py``: the Dickman rho, the
Buchstab B(u) = u*omega(u) and the Friedlander S_alpha.

They check mathematical properties: the closed forms on the first
intervals, the enclosure across precisions, continuity at the
breakpoints with the jumps each function is meant to have, the lazy
extension of the grid, and what the constructors accept. The one
comparison that is not an enclosure check, against Sage's
floating-point ``dickman_rho``, says so in its name and in its
message.

The module is loaded on its own, through ``dde_namespace()``, not through the
library.
"""


def solutions():
    """The three classes, as a dict."""
    namespace = dde_namespace()
    return {name: namespace[name]
            for name in ("DickmanRho", "BuchstabB", "FriedlanderS")}


def build(name, precision, alpha=None):
    """An instance of ``name``; ``alpha`` is required by FriedlanderS alone."""
    cls = solutions()[name]
    return cls(precision=precision) if alpha is None else cls(alpha, precision=precision)


# ------------------------------------------------------------ closed forms
#
# rho = 1 on [0, 1] and 1 - log(u) on [1, 2]; B = 0 below 1, 1 on [1, 2] --
# omega being 1/u there -- and 1 + log(u-1) on [2, 3]. The closed form is
# evaluated fifty bits higher, and has to be enclosed.

DICKMAN_CLOSED_FORM = [
    (QQ(1)/100, lambda R, x: R(1)), (QQ(1)/2, lambda R, x: R(1)),
    (QQ(99)/100, lambda R, x: R(1)), (QQ(1), lambda R, x: R(1)),
    (QQ(101)/100, lambda R, x: 1 - R(x).log()), (QQ(3)/2, lambda R, x: 1 - R(x).log()),
    (QQ(199)/100, lambda R, x: 1 - R(x).log()), (QQ(2), lambda R, x: 1 - R(x).log()),
]

BUCHSTAB_CLOSED_FORM = [
    (QQ(1)/2, lambda R, x: R(0)), (QQ(99)/100, lambda R, x: R(0)),
    (QQ(1), lambda R, x: R(1)), (QQ(101)/100, lambda R, x: R(1)),
    (QQ(3)/2, lambda R, x: R(1)), (QQ(199)/100, lambda R, x: R(1)),
    (QQ(2), lambda R, x: 1 + R(x - 1).log()),
    (QQ(201)/100, lambda R, x: 1 + R(x - 1).log()),
    (QQ(5)/2, lambda R, x: 1 + R(x - 1).log()),
    (QQ(299)/100, lambda R, x: 1 + R(x - 1).log()),
    (QQ(3), lambda R, x: 1 + R(x - 1).log()),
]


def closed_form_test(name, points, precision):
    def test():
        solution = build(name, precision)
        exact_field = RealBallField(precision + 50)
        for point, closed_form in points:
            value = solution(point)
            expected = closed_form(exact_field, point)
            assert value.overlaps(expected), (
                f"{name}(precision={precision})({point}) = {value} does not "
                f"contain the closed form {expected}")
    return test


# ------------------------------------------------- enclosure across precisions
#
# Two enclosures of the same number have to meet.

ENCLOSURE_POINTS = (QQ(5)/2, QQ(9)/2, QQ(13)/2, QQ(17)/2, QQ(23)/2)


def enclosure_test(name, alpha=None, points=ENCLOSURE_POINTS):
    def test():
        low = build(name, 30, alpha)
        high = build(name, 200, alpha)
        for point in points:
            coarse, fine = low(point), high(point)
            assert coarse.overlaps(fine), (
                f"{name} at 30 bits gives {coarse} at u = {point}, which does "
                f"not meet {fine} at 200 bits")
    return test


# ------------------------------------------------ continuity at the breakpoints
#
# Consecutive approximants have to agree at the step between them, up to the
# jump the function is meant to have there: +1 at 1 for B, +1 at alpha and -1
# at 1 for S_alpha.

def continuity_test(name, alpha=None, jumps=None, reach=8):
    jumps = jumps or {}

    def test():
        solution = build(name, 53, alpha)
        solution(reach)
        pieces = solution.approximants_list()
        for (left, step, before), (_, _, after) in zip(pieces, pieces[1:]):
            expected = before(1) + jumps.get(step, 0)
            assert after(-1).overlaps(expected), (
                f"{name} jumps at u = {step}: the approximant on [{left}, "
                f"{step}] ends at {before(1)}, the next starts at {after(-1)}, "
                f"where {expected} was expected")
    return test


# -------------------------------------------------------- the jumps of S_alpha

def jump_test(alpha, epsilon=QQ(1)/10**12):
    def test():
        solution = build("FriedlanderS", 80, alpha)
        at_alpha = solution(alpha + epsilon) - solution(alpha - epsilon)
        assert at_alpha.overlaps(RealBallField(80)(1)), (
            f"S_{alpha} jumps by {at_alpha} at alpha, where 1 was expected")
        # Away from the jump S_alpha is smooth, so the difference over the
        # window of width 2*epsilon differs from -1 by O(epsilon).
        at_one = solution(1 + epsilon) - solution(1 - epsilon)
        error = (at_one + 1).abs().upper()
        assert error < 10 * epsilon, (
            f"S_{alpha} jumps by {at_one} at 1, which is {error} away from -1, "
            f"more than the {10 * epsilon} the derivative accounts for")
    return test


# ------------------------------------------------------------ lazy extension
#
# Evaluating beyond the last step extends the grid, filling in the new
# intervals only. What was already computed is not recomputed, so an extended
# instance holds exactly the coefficients of a fresh one.

def extension_test(name, alpha=None, reach=20):
    def test():
        extended = build(name, 53, alpha)
        before = len(extended.approximants_list())
        extended(reach // 2)
        value = extended(reach)
        fresh = build(name, 53, alpha)
        expected = fresh(reach)
        assert len(extended.approximants_list()) > before, (
            f"{name}({reach}) did not extend the grid")
        assert value.overlaps(expected), (
            f"{name} extended to {reach} gives {value}, a fresh instance "
            f"{expected}")
        for (start, _, mine), (_, _, theirs) in zip(extended.approximants_list(),
                                                    fresh.approximants_list()):
            assert mine.degree() == theirs.degree(), (
                f"{name} recomputed the interval at {start}: degree "
                f"{mine.degree()} extended against {theirs.degree()} fresh")
            for degree, (a, b) in enumerate(zip(mine.list(), theirs.list())):
                assert a.mid() == b.mid() and a.rad() == b.rad(), (
                    f"{name} recomputed the interval at {start}: its "
                    f"coefficient of degree {degree} is {a} extended and {b} "
                    f"fresh")
    return test


# ----------------------------------------------- sanity against Sage's rho
#
# dickman_rho is a floating-point computation, with no error bound. We
# check that the difference is at most 1e-12 at 53 bits of requested
# precision.

SAGE_RHO_POINTS = (2, 3, QQ(9)/2, 7, 10)


def sage_dickman_sanity():
    solution = build("DickmanRho", 53)
    for point in SAGE_RHO_POINTS:
        ours = solution(point)
        theirs = RealBallField(53)(dickman_rho(RealField(53)(point)))
        difference = (ours - theirs).abs().upper()
        assert difference < 1e-12, (
            f"rho({point}) = {ours} is {difference} away from Sage's "
            f"floating-point {theirs}; this is a sanity check, not an "
            f"enclosure: neither value bounds the other")


# --------------------------------------------------------- what is accepted

def alpha_rejected_test():
    FriedlanderS = solutions()["FriedlanderS"]
    # QQ(None) is 0
    for alpha in (QQ(3)/2, QQ(1), QQ(0), -QQ(1)/2, None):
        try:
            FriedlanderS(alpha, precision=20)
        except ValueError:
            continue
        raise AssertionError(
            f"FriedlanderS({alpha}) was built: alpha has to satisfy "
            f"0 < alpha < 1, or the interval at 1 has no predecessor at "
            f"distance alpha and is left without an approximant")
    for alpha in (sqrt(2), "x", I):
        try:
            FriedlanderS(alpha, precision=20)
        except TypeError:
            continue
        raise AssertionError(f"FriedlanderS({alpha}) was built: alpha has to "
                             f"be rational")


def arguments_accepted_test():
    for name, alpha in (("DickmanRho", None), ("BuchstabB", None),
                        ("FriedlanderS", QQ(2)/5)):
        solution = build(name, 30, alpha)
        for point, others in (
                (QQ(5)/2, (RR(2.5), float(2.5), RealBallField(53)(2.5))),
                (QQ(3), (3, int(3), RR(3), RealBallField(53)(3)))):
            reference = solution(point)
            for argument in others:
                value = solution(argument)
                assert value.overlaps(reference), (
                    f"{name} at u = {point} given as "
                    f"{type(argument).__name__} gives {value}, against "
                    f"{reference} given as a Rational")


# ------------------------------------------------------------------ registry

TESTS = {
    "dde/closedform/dickman": (
        closed_form_test("DickmanRho", DICKMAN_CLOSED_FORM, 53), ("quick",)),
    "dde/closedform/dickman-150": (
        closed_form_test("DickmanRho", DICKMAN_CLOSED_FORM, 150), ("quick",)),
    "dde/closedform/buchstab": (
        closed_form_test("BuchstabB", BUCHSTAB_CLOSED_FORM, 53), ("quick",)),
    "dde/closedform/buchstab-150": (
        closed_form_test("BuchstabB", BUCHSTAB_CLOSED_FORM, 150), ("quick",)),
    "dde/enclosure/dickman": (enclosure_test("DickmanRho"), ("quick",)),
    "dde/enclosure/buchstab": (enclosure_test("BuchstabB"), ("quick",)),
    "dde/enclosure/friedlander": (
        enclosure_test("FriedlanderS", QQ(2)/5,
                       (QQ(3)/4,) + ENCLOSURE_POINTS), ("quick",)),
    "dde/continuity/dickman": (continuity_test("DickmanRho"), ("quick",)),
    "dde/continuity/buchstab": (
        continuity_test("BuchstabB", jumps={QQ(1): 1}), ("quick",)),
    "dde/continuity/friedlander": (
        continuity_test("FriedlanderS", QQ(2)/5,
                        jumps={QQ(2)/5: 1, QQ(1): -1}), ("quick",)),
    "dde/jumps/friedlander-2over5": (jump_test(QQ(2)/5), ("quick",)),
    "dde/jumps/friedlander-1over2": (jump_test(QQ(1)/2), ("quick",)),
    "dde/extension/dickman": (extension_test("DickmanRho"), ("quick",)),
    "dde/extension/buchstab": (extension_test("BuchstabB"), ("quick",)),
    "dde/extension/friedlander": (
        extension_test("FriedlanderS", QQ(2)/5), ("quick",)),
    "dde/sanity/dickman-against-sage": (sage_dickman_sanity, ("quick",)),
    "dde/errors/friedlander-alpha": (alpha_rejected_test, ("quick",)),
    "dde/arguments/accepted": (arguments_accepted_test, ("quick",)),
}
