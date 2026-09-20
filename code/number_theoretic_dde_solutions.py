r"""
Three particular solutions of differential-delay equations of interest in number theory

Each of them is represented as a piecewise polynomial whose coefficients lie in
a ``RealBallField``, so that every value returned is a rigorous enclosure:

- :class:`DickmanRho`, the Dickman `\rho` function;
- :class:`BuchstabB`, the map `B(u) = u \omega(u)` with `\omega` the Buchstab
  function;
- :class:`FriedlanderS`, a variant of the Friedlander `\sigma` function.

The three share the base class :class:`DifferentialDelaySolution`, which holds
the Marsaglia-Zaman-Marsaglia scheme they all follow; a subclass describes the
differential equation it solves, its known first intervals, and, if they are not
the non-negative integers, its steps.

This file can be ``load``ed in a Sage session, as ``sieve_integral.py`` does, or
imported as a module. Beware that a ``load``ed file does not go through the Sage
preparser, so that ``1 / 2`` in it would be a Python float: the grid of steps is
built in ``QQ`` for that reason.
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


# Everything is taken from sage.all: importing a Sage submodule, such as
# sage.functions.other for ceil, before sage.all in a fresh Python process
# trips a circular import inside Sage itself.
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

    """

    def __init__(self, precision=53):
        r"""
        INPUT:

        - ``precision`` (default: 53) -- the precision, in bits, of the
        ``RealBallField`` the approximants are defined over.
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
    # out by the functional equations
    # _scalar_ring : RealBallField
    # _polynomial_ring : Polynomial Ring over _scalar_ring of 1 variable
    # _approximants : dictionary of elements of _polynomial_ring. The element with key s0 is a
    # polynomial which approximates the map F
    #      over the interval starting with s0, mapped to [-1, 1], so
    # that polynomial(-1) approximates F(step_begin).
    # _taylor_degree : number representing the largest degree we
    # perform Taylor approximation at.

    # ------------------------------------------------------------------
    # The grid of steps
    # ------------------------------------------------------------------

    def _compute_steps(self, smax):
        """
        Computes the list _steps of the points delimiting the intervals where F
        is smooth. By default these are the integers up to smax; a subclass
        whose grid is different overrides this.

        The steps are rationals, and not Python integers, because this module is
        also loaded as a plain .py file, which the Sage preparser does not see:
        (a + b) / 2 on two Python integers would return a float.
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

        The last step begins no interval, and the seeds are in _approximants
        already, so this is exactly what _compute_approximants has to run over.
        An approximant is never recomputed: recomputing one returns the very
        same coefficients, midpoint and radius, since it is built from the
        stored -- already truncated -- approximant of the interval before it.
        Extending the grid therefore only fills in the new intervals.
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

        A vanishing denominator_minorant would make the degree below infinite;
        the callers which have one, at the left end of the interval where the
        delayed term is still zero, are covered by the test on the numerator.
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
        function, which drops the same terms (truncate_abs) but, being a
        floating-point computation, does not account for them.

        The discarded coefficients are accumulated as a *radius*, through
        add_error, and not as a value: the tail they stand for is only known to
        be at most their sum in absolute value, so adding that sum would move
        the approximant where it has to be widened, and the resulting ball
        could then fail to contain the true value.
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

    AUTHORS:

    - Sary Drappeau (2026-06)

    AI DISCLOSURE:

    GPT 5.2 provided the basic structure of the class.

    REFERENCES:

    - G. Marsaglia, A. Zaman, J. Marsaglia. "Numerical
      Solutions to some Classical Differential-Difference Equations."
      Mathematics of Computation, Vol. 53, No. 187 (1989).

    """

    def __init__(self, alpha, precision=53):
        # alpha has to be rational: the steps m + n * alpha are compared for
        # equality as the grid is built, so they have to be exact. A real
        # argument is accepted, but QQ reconstructs a nearby rational from it
        # by continued fractions, which gives 2/5 for 0.4 but a denominator in
        # the millions for a number which is not close to a simple fraction,
        # and as many steps; pass a Rational to know what is computed.
        try:
            alpha = QQ(alpha)
        except (TypeError, ValueError) as exc:
            raise TypeError("the parameter alpha should be rational") from exc
        if not 0 < alpha < 1:
            # For alpha >= 1, the interval starting at 1 has no predecessor at
            # distance alpha, so the recursion never reaches it and its
            # approximant would be left undefined.
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
                # Below u = 1 there is no delayed term at distance 1 at all.
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
