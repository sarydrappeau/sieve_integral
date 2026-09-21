r"""
Three particular solutions of differential-delay equations of interest in number theory

Each of them is represented as a piecewise polynomial whose coefficients lie in
a ``RealBallField`` (rigorous enclosure)

- :class:`DickmanRho`, the Dickman `\rho` function;
- :class:`BuchstabB`, the map `B(u) = u \omega(u)` with `\omega` the Buchstab
  function;
- :class:`FriedlanderS`, a variant of the Friedlander `\sigma` function.

The three share the base class :class:`DifferentialDelaySolution`,
written with the same Marsaglia-Zaman-Marsaglia scheme ; a subclass
describes the differential equation it solves, its known first
intervals, and, if they are not the non-negative integers, its steps.

EXAMPLES::

    sage: DickmanRho(precision=53)(2).overlaps(RBF(1 - log(2)))
    True
    sage: BuchstabB(precision=53)(3).overlaps(RBF(1 + log(2)))
    True
    sage: FriedlanderS(2/5, precision=53)(1/2)
    1.000000000000000

AUTHORS:

- Sary Drappeau (2026-06)

AI DISCLOSURE:

ChatGPT 5.2 provided the basic structure of the class.
"""
# ****************************************************************************
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  https://www.gnu.org/licenses/
# ****************************************************************************


from sage.all import QQ, RR, PolynomialRing, RealBallField, ceil, floor


class DifferentialDelaySolution:
    r"""
    Base class for the piecewise real-analytic solutions of differential-delay
    equations implemented in this module.

    Such a solution `F` is smooth away from a discrete set of *steps*, and on
    the interval between two consecutive steps it is represented by a polynomial
    with coefficients in a ``RealBallField``. The polynomials are normalized so
    that their variable ranges over `[-1, 1]`: if `P` is the approximant
    attached to the interval `[s_i, s_{i+1}]`, then `P(-1)` approximates
    `F(s_i)` and `P(1)` approximates `F(s_{i+1})`.

    The computation follows the Marsaglia-Zaman-Marsaglia method, in the form
    used by R. Bradshaw's ``dickman_rho`` implementation in Sage, the difference
    being that the coefficients live in a ``RealBallField``, so that the
    approximation is rigorous.

    A subclass provides:

    - ``_seed_approximants()``, the approximants known outright, on the first
      intervals, where the delay term vanishes;

    - ``_compute_approximants()``, which integrates the differential equation
      from one interval to the next;

    - ``_compute_steps()``, only if the steps are not the non-negative integers.

    REFERENCES:

    - G. Marsaglia, A. Zaman, J. Marsaglia. "Numerical
      Solutions to some Classical Differential-Difference Equations."
      Mathematics of Computation, Vol. 53, No. 187 (1989).

    TESTS:

    Two enclosures of the same number have to meet, so the value returned at a
    given precision must overlap the value returned at a higher one. This is
    what fails when the coefficients discarded by ``_truncate`` are accumulated
    as a value instead of as a radius, which moves the ball instead of widening
    it::

        sage: low, high = BuchstabB(precision=30), BuchstabB(precision=200)
        sage: all(low(u).overlaps(high(u)) for u in (5/2, 13/2, 15/2, 21/2))
        True
        sage: low, high = DickmanRho(precision=30), DickmanRho(precision=200)
        sage: all(low(u).overlaps(high(u)) for u in (5/2, 11/2, 15/2))
        True
        sage: low = FriedlanderS(2/5, precision=30)
        sage: high = FriedlanderS(2/5, precision=200)
        sage: all(low(u).overlaps(high(u)) for u in (3/4, 11/5, 17/5, 4))
        True

    """

    def __init__(self, precision=53):
        r"""
        INPUT:

        - ``precision`` (default: 53) -- the precision, in bits, of the
        ``RealBallField`` the approximants are defined over.

        EXAMPLES::

            sage: DickmanRho(precision=30).base_ring()
            Real ball field with 30 bits of precision
        """
        self._taylor_degree = floor(precision * RR(2).log() / RR(3).log())
        self._scalar_ring = RealBallField(precision)
        self._polynomial_ring = PolynomialRing(self._scalar_ring, "t")

        # Compute a few values
        self._compute_steps(3)

        self._approximants = dict(self._seed_approximants())
        self._compute_approximants()

    # local variables :
    # _steps : list of rational numbers, delimiting the intervals cut
    #      out by the functional equations
    # _scalar_ring : RealBallField
    # _polynomial_ring : Polynomial Ring over _scalar_ring of 1 variable
    # _approximants : dictionary of elements of _polynomial_ring. The element with key s0 is a
    #      polynomial which approximates the map F
    #      over the interval starting with s0, mapped to [-1, 1], so
    #      that polynomial(-1) approximates F(step_begin).
    # _taylor_degree : number representing the largest degree we
    #      perform Taylor approximation at.

    # ------------------------------------------------------------------
    # The grid of steps
    # ------------------------------------------------------------------

    def _compute_steps(self, smax):
        """
        Computes the list _steps of the points delimiting the intervals where F
        is smooth. By default these are the integers up to smax. A subclass
        whose grid is different overrides this.
        """
        self._steps = [QQ(i) for i in range(floor(smax) + 1)]

    def _successor(self, step):
        """
        Given an element step of _steps, returns the element which
        immediately follows step if it exists, otherwise returns
        None.
        """
        if step not in self._steps or step == self._steps[-1]:
            return None
        return self._steps[self._steps.index(step) + 1]

    def _centre(self, current_left):
        """
        Returns the center of the interval for which current_left is the left end.
        """
        current_right = self._successor(current_left)
        if current_right is None:
            return None
        return (current_left + current_right) / 2

    def _predecessor_leq(self, target):
        """
        Returns the largest step which is less than or equal to target, or None
        if there is none.
        """
        candidates = [step for step in self._steps if step <= target]
        return max(candidates) if candidates else None

    def _to_unit_interval(self, current_left, x):
        """
        Maps the interval starting with current_left onto [-1, 1].
        """
        current_center = self._centre(current_left)
        current_radius = current_center - current_left
        return (x - current_center) / current_radius

    def _from_unit_interval(self, current_left, x):
        """
        Inverse of _to_unit_interval.
        """
        current_center = self._centre(current_left)
        current_radius = current_center - current_left
        return current_center + x * current_radius

    # ------------------------------------------------------------------
    # The approximants
    # ------------------------------------------------------------------

    def _seed_approximants(self):
        """
        Returns a dictionary of the approximants known outright, on the first
        intervals. To be provided by the subclass.
        """
        raise NotImplementedError

    def _compute_approximants(self):
        """
        Fills in _approximants by integrating the differential equation from one
        interval to the next, over the steps _steps_to_compute() returns. To be
        provided by the subclass.
        """
        raise NotImplementedError

    def _steps_to_compute(self):
        """
        The steps whose interval has no approximant yet, in increasing order.

        EXAMPLES::

            sage: rho = DickmanRho(precision=30)
            sage: rho._steps_to_compute()
            []
            sage: rho._compute_steps(6)
            sage: rho._steps_to_compute()
            [3, 4, 5]
        """
        return [step for step in self._steps[:-1]
                if step not in self._approximants]

    def _adjust_constant(self, current_left):
        """
        Add a number to the polynomial approximation on the interval
        starting with current_left, so that its value coincides with the value
        of the polynomial which precedes it, evaluated at the same
        point.
        """
        i = self._steps.index(current_left)
        if i == 0:
            return

        previous_left = self._steps[i - 1]

        # recall that the polynomials are normalized so that their
        # variables vary over [-1, 1]
        val_left = self._approximants[previous_left](1)
        val_right = self._approximants[current_left](-1)

        correction = val_left - val_right
        self._approximants[current_left] += correction

    def _expand_quotient(self, numerator, denominator, denominator_minorant):
        """
        Compute a Taylor approximation at 0, on the interval [-1, 1],
        of numerator / denominator, where numerator, denominator are
        polynomials with denominator linear (deg denominator = 1). It
        is assumed that denominator_minorant is a lower-bound for
        |denominator|. The method is borrowed from R. Bradshaw's
        sagemath implementation of the Dickman rho function.
        """
        if numerator == self._polynomial_ring(0):
            return self._polynomial_ring(0)
        if denominator.degree() != 1:
            raise NotImplementedError
        degree = ceil(max(RR(0), -denominator_minorant.log().lower())
                      * RR(2).log() / RR(3).log()
                      + self._taylor_degree)
        # Let numerator_reversed, denominator_reversed be the reverses
        # of numerator, denominator.  A quick computation reveals that
        # numerator/denominator = quotient + x^N c /
        # denominator, where numerator_reversed / denominator_reversed
        # = quotient_reversed + c is the Euclidean division, and
        # quotient_reversedr is the reversed of quotient_reversed.
        numerator_reversed = numerator.reverse(degree)
        denominator_reversed = denominator.reverse()
        quotient_reversed, remainder = numerator_reversed.quo_rem(denominator_reversed)
        quotient = quotient_reversed.reverse()
        remainder = remainder.coefficient(0)
        return quotient + quotient.base_ring()(0).add_error((remainder / denominator_minorant).above_abs().upper())

    def _truncate(self, current_left):
        """
        Having computed the polynomial approximation P on the interval starting
        with current_left, truncates the terms of high degree whose size is
        individually below the error already carried by the constant coefficient
        of P, and adds their total to that error. Their number is not bounded, so
        the error can come out a few times larger than the one it is compared
        against; this is a deliberate trade of width against degree. It is also
        borrowed from R. Bradshaw's sagemath implementation of the Dickman rho
        function, which drops the same terms (truncate_abs).
        """
        approximant = self._approximants[current_left]
        total_error = approximant.base_ring()(0)
        constant_error = approximant.coefficient(0).rad()
        truncation_degree = approximant.degree()
        while truncation_degree >= 0 and approximant.coefficient(truncation_degree).abs().upper() < constant_error:
            total_error = total_error.add_error(approximant.coefficient(truncation_degree).above_abs())
            truncation_degree -= 1
        self._approximants[current_left] = approximant.truncate(truncation_degree + 1) + total_error

    # ------------------------------------------------------------------
    # Evaluation
    # ------------------------------------------------------------------

    def __call__(self, sval):
        """
        Returns a RealBallField element which approximates F(u) with u = sval.

        The argument may be of any real type, and one beyond the last step
        extends the grid.

        EXAMPLES::

            sage: rho = DickmanRho(precision=30)
            sage: all(rho(u).overlaps(rho(2)) for u in (2, QQ(2), RR(2), float(2), RBF(2)))
            True
            sage: rho(-1)
            0
            sage: rho(12).upper() < 1e-5
            True
            sage: rho._steps[-1]
            14
        """
        scalar_ring = self._scalar_ring

        if sval >= self._steps[-1]:
            # Take a small margin
            self._compute_steps(RR(sval) * 1.1 + 1)
            self._compute_approximants()

        if sval <= 0:
            return scalar_ring(0)

        current_left = self._predecessor_leq(sval)
        if current_left is None:
            return scalar_ring(0)

        # sval is coerced into the ball field before the polynomial is
        # evaluated, so that an argument of any real type is accepted.
        return self._approximants[current_left](
            self._to_unit_interval(current_left, scalar_ring(sval)))

    def base_ring(self):
        """
        Returns the RealBallField the approximants are defined over.

        EXAMPLES::

            sage: BuchstabB(precision=20).base_ring()
            Real ball field with 20 bits of precision
        """
        return self._scalar_ring

    def approximants_list(self):
        """
        This is the function which is used most for the sieve
        integrals procedures. Returns a list of tuples (start, end,
        polynomial), where polynomial is a polynomial with
        coefficients in a RealBallField which approximates F on the
        interval [start, end] mapped to [-1, 1], so that polynomial(1)
        = F(end) for instance.

        EXAMPLES::

            sage: B = BuchstabB(precision=53)
            sage: [(start, end) for start, end, approximant in B.approximants_list()]
            [(0, 1), (1, 2), (2, 3)]
            sage: start, end, approximant = B.approximants_list()[2]
            sage: approximant(1).overlaps(RBF(1 + log(2)))
            True
        """
        return [[self._steps[i], self._steps[i + 1], self._approximants[self._steps[i]]] for i in range(len(self._steps) - 1)]


class FriedlanderS(DifferentialDelaySolution):
    r"""
    The Friedlander `S_\alpha` function is the unique function on
    `\mathbb R` which is continuous on `D =
    \mathbb{R}\setminus\{\alpha, 1\}`, real-analytic away
    from numbers of the shape `m + n\alpha` (`m,
    n\in\mathbb{Z}_{\geq 0}`), with a jump of size `1` at
    `\alpha` and a jump of size `-1` at `1`, satisfying the
    differential equation

    .. MATH::

        S_\alpha'(x) = \frac{S_\alpha(x-\alpha)}{x-\alpha}
        - \frac{S_\alpha(x-1)}{x-1}, \qquad (x\neq \alpha, 1).

    The map `\sigma(u, v) = 1/v S_{u/v}(u)` corresponds to
    Friedlander's sigma function [Friedlander1976].

    INPUT:

    - ``alpha`` -- a rational number with `0 < \alpha < 1`.

    - ``precision`` (default: 53) -- the precision, in bits, of the
    ``RealBallField`` the approximants are defined over.

    EXAMPLES::

        sage: S = FriedlanderS(1/2, precision=53)
        sage: S(1/4)
        0
        sage: S(3/4)
        1.000000000000000

    The jump at `\alpha` has size `1`, and the one at `1` has size `-1`::

        sage: S = FriedlanderS(2/5, precision=80)
        sage: eps = QQ(1)/10^12
        sage: (S(2/5 + eps) - S(2/5 - eps) - 1).abs().upper() < 1e-20
        True
        sage: (S(1 + eps) - S(1 - eps) + 1).abs().upper() < 1e-9
        True

    Friedlander proved that `\sigma(u, v)` tends to `e^{-\gamma}\rho(u)` as
    `v \to \infty`, with `\rho` the Dickman function. At `u = 3.1`, over the
    range of `v` the enclosure can still resolve, the relative gap is at most
    `7/v`::

        sage: u = QQ(31)/10
        sage: target = RBF(-euler_gamma).exp() * DickmanRho(precision=53)(u)
        sage: gaps = [(v, 1 - FriedlanderS(u/v, precision=53)(u) / (v * target))
        ....:         for v in (20, 40, 80, 160)]
        sage: all((gap.abs() * v).upper() < 7 for v, gap in gaps)
        True

    The parameter has to be rational, with `0 < \alpha < 1`; a real one is
    reconstructed as a nearby rational, by continued fractions::

        sage: FriedlanderS(3/2)
        Traceback (most recent call last):
        ...
        ValueError: the parameter alpha should satisfy 0 < alpha < 1
        sage: FriedlanderS(sqrt(2))
        Traceback (most recent call last):
        ...
        TypeError: the parameter alpha should be rational
        sage: FriedlanderS(0.4).alpha
        2/5

    """

    def __init__(self, alpha, precision=53):
        try:
            alpha = QQ(alpha)
        except (TypeError, ValueError) as exc:
            raise TypeError("the parameter alpha should be rational") from exc
        if not 0 < alpha < 1:
            raise ValueError("the parameter alpha should satisfy "
                             "0 < alpha < 1")
        self.alpha = alpha
        super().__init__(precision)

    # local variables :
    # alpha : rational, value of the parameter (cf. description above)

    def _compute_steps(self, smax):
        """
        Computes the list _steps of all numbers of the shape m + n * alpha at most smax.
        These values delimit the various intervals where F is smooth.
        """
        steps = set()
        alpha = self.alpha
        i = 0
        while i <= smax:
            j = 0
            while i + j * alpha <= smax:
                steps.add(i + j * alpha)
                j += 1
            i += 1
        self._steps = sorted(steps)

    def _seed_approximants(self):
        # S vanishes below alpha, and the jump of size 1 at alpha leaves it
        # equal to 1 on the interval which starts there.
        return {QQ(0): self._polynomial_ring(0),
                self.alpha: self._polynomial_ring(1)}

    def _compute_approximants(self):
        r"""
        Integrates `S_\alpha'(u) = S_\alpha(u-\alpha)/(u-\alpha) -
        S_\alpha(u-1)/(u-1)` from one interval to the next.
        """
        alpha = self.alpha
        scalar_ring = self._scalar_ring
        poly_gen = self._polynomial_ring.gen()

        for current_left in self._steps_to_compute():
            current_radius = self._centre(current_left) - current_left
            # u, as a function of the variable of the approximant to be
            # computed, which ranges over [-1, 1]
            u = self._from_unit_interval(current_left, poly_gen)

            delayed_alpha = self._predecessor_leq(current_left - alpha)
            term = self._expand_quotient(
                self._approximants[delayed_alpha](
                    self._to_unit_interval(delayed_alpha, u - alpha)),
                u - alpha,
                scalar_ring(current_left - alpha))

            delayed_one = self._predecessor_leq(current_left - 1)
            if delayed_one is not None:
                # Below u = 1 there is no delayed term at all.
                term -= self._expand_quotient(
                    self._approximants[delayed_one](
                        self._to_unit_interval(delayed_one, u - 1)),
                    u - 1,
                    scalar_ring(current_left - 1))

            self._approximants[current_left] = current_radius * term.integral()
            self._adjust_constant(current_left)
            if current_left == 1:
                # Account for the discontinuity at u = 1.
                self._approximants[current_left] -= 1
            self._truncate(current_left)


class DickmanRho(DifferentialDelaySolution):
    r"""
    The Dickman `\rho` function, the continuous solution of

    .. MATH::

        u \rho'(u) = -\rho(u-1), \qquad (u > 1),

    which is equal to `1` on `[0, 1]`.

    Performs roughly the same operations as R. Bradshaw's
    implementation dickman_rho, the main difference being that it
    returns a RealBallField element, so that the approximation is
    rigorous.

    EXAMPLES::

        sage: rho = DickmanRho(precision=53)
        sage: rho(1/2)
        1.000000000000000
        sage: rho(-1)
        0

    On `[1, 2]` one has `\rho(u) = 1 - \log u`::

        sage: rho(2).overlaps(RBF(1 - log(2)))
        True
        sage: rho(3/2).overlaps(RBF(1 - log(3/2)))
        True

    At `u = 30`, where `\rho(u) \approx 3.27 \cdot 10^{-50}`, the enclosure
    contains the value Sage's own floating-point ``dickman_rho`` gives, and the
    two agree to some thirty digits, which is the width of the enclosure at 300
    bits. Sage's value is not rigorous, so this is a cross-check and not an
    enclosure test::

        sage: ours = DickmanRho(precision=300)(30)
        sage: theirs = RealBallField(300)(dickman_rho(RealField(300)(30)))
        sage: ours.overlaps(theirs)
        True
        sage: ((ours - theirs) / theirs).abs().upper() < 1e-30
        True

    """

    def _seed_approximants(self):
        # rho = 1 on [0, 1].
        return {QQ(0): self._polynomial_ring(1)}

    def _compute_approximants(self):
        r"""
        Integrates `\rho'(u) = -\rho(u-1)/u` from one interval to the next.
        """
        scalar_ring = self._scalar_ring
        poly_gen = self._polynomial_ring.gen()

        for current_left in self._steps_to_compute():
            previous_left = current_left - 1

            # On the interval starting with current_left, u = centre + x/2 with
            # x in [-1, 1], and u - 1 is described by that same x on the
            # interval before. The denominator u is at least current_left there.
            term = self._expand_quotient(
                self._approximants[previous_left](poly_gen),
                poly_gen / 2 + scalar_ring(self._centre(current_left)),
                scalar_ring(current_left))
            self._approximants[current_left] = (
                scalar_ring(-QQ(1) / 2) * term.integral())
            self._adjust_constant(current_left)
            self._truncate(current_left)


class BuchstabB(DifferentialDelaySolution):
    r"""
    An element in this class represents the map  B(u) = u ω(u)  where
    ω is the classical Buchstab function.

    It is the continuous solution of

    .. MATH::

        B'(u) = \frac{B(u-1)}{u-1}, \qquad (u > 2),

    which is equal to `1` on `[1, 2]`, and to `0` below `1`.

    EXAMPLES::

        sage: B = BuchstabB(precision=53)
        sage: B(1/2)
        0
        sage: B(3/2)
        1.000000000000000

    On `[2, 3]` one has `B(u) = 1 + \log(u-1)`::

        sage: B(3).overlaps(RBF(1 + log(2)))
        True
        sage: B(5/2).overlaps(RBF(1 + log(3/2)))
        True

    As `u \to \infty`, `\omega(u) = B(u)/u` tends to `e^{-\gamma}`, and at
    `u = 30` the two already agree to well below the working precision::

        sage: target = RealBallField(200)(-euler_gamma).exp()
        sage: omega30 = BuchstabB(precision=200)(30) / 30
        sage: omega30.overlaps(target)
        True
        sage: (omega30 - target).abs().upper() < 1e-50
        True

    """

    def _seed_approximants(self):
        # B vanishes below 1, and B = 1 on [1, 2].
        return {QQ(0): self._polynomial_ring(0),
                QQ(1): self._polynomial_ring(1)}

    def _compute_approximants(self):
        r"""
        Integrates `B'(u) = B(u-1)/(u-1)` from one interval to the next.
        """
        scalar_ring = self._scalar_ring
        poly_gen = self._polynomial_ring.gen()

        for current_left in self._steps_to_compute():
            previous_left = current_left - 1

            # On the interval starting with current_left, u - 1 is described by
            # the same variable x in [-1, 1] as u, on the interval before: it is
            # the centre of that one plus x/2, and at least previous_left.
            term = self._expand_quotient(
                self._approximants[previous_left](poly_gen),
                poly_gen / 2 + scalar_ring(self._centre(previous_left)),
                scalar_ring(previous_left))
            self._approximants[current_left] = (
                scalar_ring(QQ(1) / 2) * term.integral())
            self._adjust_constant(current_left)
            self._truncate(current_left)
