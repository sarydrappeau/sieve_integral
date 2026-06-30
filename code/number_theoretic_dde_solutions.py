"""
Three particular solutions of differential-delay equations of interest in number theory
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


from sage.all import RR, PolynomialRing, floor
from sage.rings.real_arb import RealBallField
from sage.functions.other import ceil


class FriedlanderS:
    """
    The Friedlander `S_\\alpha` function is the unique fonction on
    `\\mathbb R` which is continuous on `D =
    \\mathbb{R}\\setminus\\{\\alpha, 1\\}`, real-analytic away
    from numbers of the shape `m + n\\alpha` (`m,
    n\\in\\mathbb{Z}_{\\geq 0}`, with a jump of size `1` at
    `\\alpha` and a jump of size `-1` at `1` satisfying the
    differential equation

    .. MATH::

    S_\\alpha'(x) = \\frac{S_\\alpha(x-\\alpha)}{x-\\alpha)
    - \\frac{S_alpha(x-1)}{x-1}, \\qquad (x\\neq \\alpha, 1).

    The map `\\sigma(u, v) = 1/v S_{u/v}(u)` corresponds to
    Friedlander's sigma function [Friedlander1976].

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
        self.alpha = alpha
        self._taylor_degree = floor(precision * RR(2).log() / RR(3).log())

        scalar_ring = RealBallField(precision)
        polynomial_ring = PolynomialRing(scalar_ring, "t")
        self._scalar_ring = scalar_ring
        self._polynomial_ring = polynomial_ring

        # Compute a few values
        self._compute_steps(3)

        approximants = {}
        approximants[0] = polynomial_ring(0)
        approximants[alpha] = polynomial_ring(1)
        self._approximants = approximants

        self._compute_approximants()

    # local variables :
    # alpha : rational, value of the parameter (cf. description above)
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

    def _compute_steps(self, smax):
        """
        Computes the list _steps of all numbers of the shape m + n * \\alpha at most smax.
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

    def _successeur(self, step):
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
        current_right = self._successeur(current_left)
        if current_right is None:
            return None
        return (current_left + current_right) / 2

    def _predecessor_leq(self, vals, target):
        """
        Returns the largest value of the iterable vals which is less
        than or equal to target.
        """
        cand = [v for v in vals if v <= target]
        return max(cand) if cand else None

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

    def _developpe(self, numerator, denominator, denominator_minorant):
        """
        Compute a Taylor approximation at 0, on the interval [-1, 1],
        of numerator / denominator, where numerator, denominator are
        polynomials with denominator linear (deg denominator = 1). It
        is assumed that denominator_minorant is a lower-bound for
        |denominator|. The method is borrowed from
        polynomial_ring. Bradshaw's sagemath implementation of the
        Dickman rho function.
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
        Having computed the polynomial approximation P on the interval
        starting with current_left, truncates all terms of high degree so long
        as their cumulated contribution is at most the error in the
        constant coefficient of P. This is also borrowed from
        R. Bradshaw's sagemath implementation of the Dickman rho
        function.
        """
        current_approximant = self._approximants[current_left]
        total = current_approximant.base_ring()(0)
        erreur = current_approximant.coefficient(0).rad()
        truncation_degree = current_approximant.degree()
        while truncation_degree >= 0 and current_approximant.coefficient(truncation_degree).abs().upper() < erreur:
            total += current_approximant.coefficient(truncation_degree).above_abs()
            truncation_degree -= 1
        self._approximants[current_left] = current_approximant.truncate(truncation_degree + 1) + current_approximant.base_ring()(0).add_error(total)

    def _to_unit_interval(self, current_left, x):
        current_center = self._centre(current_left)
        current_radius = current_center - current_left
        return (x - current_center) / current_radius

    def _from_unit_interval(self, current_left, x):
        current_center = self._centre(current_left)
        current_radius = current_center - current_left
        return current_center + x * current_radius


    def _compute_approximants(self):
        """
        Computes the polynomial approximations. This is run once at
        the initialization, and also possibly every time a new value
        is requested for computation which lies beyond the last
        element of steps.
        """
        alpha = self.alpha
        predecessor_leq = self._predecessor_leq
        centre = self._centre
        developpe = self._developpe

        for current_left in self._steps:
            if current_left in (0, alpha, max(self._steps)):
                # On the first two intervals, the map F is constant
                continue

            minusalpha_left = predecessor_leq(self._steps, current_left - alpha)
            minusone_left = predecessor_leq(self._steps, current_left - 1)
            current_radius = centre(current_left) - current_left

            polynomial_ring = self._polynomial_ring
            scalar_ring = self._scalar_ring
            poly_gen = polynomial_ring.gen()

            if minusalpha_left in self._approximants:
                term1 = developpe(self._approximants[minusalpha_left](self._to_unit_interval(minusalpha_left,
                                                                                             self._from_unit_interval(current_left, poly_gen) - alpha)),
                                  self._from_unit_interval(current_left, poly_gen) - alpha,
                                  scalar_ring(current_left - alpha))

                if minusone_left in self._approximants:
                    term2 = developpe(self._approximants[minusone_left](self._to_unit_interval(minusone_left,
                                                                                               self._from_unit_interval(current_left, poly_gen) - 1)),
                                      self._from_unit_interval(current_left, poly_gen) - 1,
                                      scalar_ring(current_left - 1))

                    self._approximants[current_left] = current_radius * (term1 - term2).integral()
                else:
                    self._approximants[current_left] = current_radius * term1.integral()

                self._adjust_constant(current_left)
                if current_left == 1:
                    # Account for the discontinuity at u = 1.
                    self._approximants[current_left] -= 1
                    self._truncate(current_left)

    def __call__(self, sval):
        """
        Returns a RealBallField element which approximates F(u) with u = sval.
        """
        scalar_ring = self._scalar_ring

        if sval >= self._steps[-1]:
            # Take a small margin
            self._compute_steps(RR(sval) * 1.1 + 1)
            self._compute_approximants()

        steps = self._steps

        if sval <= 0:
            return scalar_ring(0)

        current_left = self._predecessor_leq(steps, sval)
        if current_left is None:
            return scalar_ring(0)

        return self._approximants[current_left](self._to_unit_interval(current_left, sval))

    def base_ring(self):
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


class DickmanRho:
    """
    Special case of FriedlanderSigma class corresponding to the
    Dickman Rho function.  Performs roughly the same operations as
    R. Bradshaw's implementation dickman_rho, the main difference
    being that it returns a RealBallField element, so that the
    approximation is rigorous.
    """

    def _compute_steps(self, smax):
        self._steps = list(range(floor(smax) + 1))

    def _successeur(self, current_left):
        if current_left not in self._steps or current_left == self._steps[-1]:
            return None
        return current_left + 1

    def _centre(self, current_left):
        return current_left + 1 / 2

    def _predecessor_leq(self, vals, target):
        cand = [v for v in vals if v <= target]
        return max(cand) if cand else None

    def _adjust_constant(self, current_left):
        if current_left in (0, len(self._steps) - 1):
            return

        previous_left = current_left - 1

        value_left = self._approximants[previous_left](1)
        value_right = self._approximants[current_left](-1)

        correction = value_left - value_right
        self._approximants[current_left] += correction

    def _developpe(self, numerator, denominator, denominator_minorant):
        """
        Compute a Taylor approximation at 0, on the interval [-1, 1],
        of numerator / denominator, where numerator, denominator are
        polynomials with denominator linear (deg denominator = 1). It
        is assumed that denominator_minorant is a lower-bound for
        |denominator|. The method is borrowed from
        polynomial_ring. Bradshaw's sagemath implementation of the
        Dickman rho function.
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
        approximant = self._approximants[current_left]
        total_error = approximant.base_ring()(0)
        constant_error = approximant.coefficient(0).rad()
        truncation_degree = approximant.degree()
        while truncation_degree >= 0 and approximant.coefficient(truncation_degree).abs().upper() < constant_error:
            total_error = total_error.add_error(approximant.coefficient(truncation_degree).above_abs())
            truncation_degree -= 1
        self._approximants[current_left] = approximant.truncate(truncation_degree + 1) + total_error

    def _compute_approximants(self):
        for current_left in self._steps:
            if current_left in (0, max(self._steps)):
                continue

            previous_left = current_left - 1

            current_center = current_left + 1 / 2

            polynomial_ring = self._polynomial_ring
            scalar_ring = self._scalar_ring
            poly_gen = polynomial_ring.gen()
            if previous_left in self._approximants:
                term = self._developpe(self._approximants[previous_left](poly_gen),
                                       poly_gen / 2 + scalar_ring(current_center),
                                       scalar_ring(-1 / 2 + current_center))
                self._approximants[current_left] = scalar_ring(-1 / 2) * term.integral()
                self._adjust_constant(current_left)
                self._truncate(current_left)

    def __init__(self, precision=53):
        self._taylor_degree = floor(precision * RR(2).log() / RR(3).log())

        scalar_ring = RealBallField(precision)
        polynomial_ring = PolynomialRing(scalar_ring, "t")
        self._scalar_ring = scalar_ring
        self._polynomial_ring = polynomial_ring

        self._compute_steps(3)

        approximants = {}
        approximants[0] = polynomial_ring(1)
        self._approximants = approximants

        self._compute_approximants()

    def __call__(self, sval):
        scalar_ring = self._scalar_ring
        if sval >= self._steps[-1]:
            self._compute_steps(RR(sval) * 1.1 + 1)
            self._compute_approximants()

        steps = self._steps

        if sval <= 0:
            return scalar_ring(0)

        current_left = self._predecessor_leq(steps, sval)
        if current_left is None:
            return scalar_ring(0)

        current_center = self._centre(current_left)

        return self._approximants[current_left](2 * (scalar_ring(sval) - scalar_ring(current_center)))

    def base_ring(self):
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


class BuchstabB:
    """
    An element in this class represents the map  B(u) = u ω(u)  where
    ω is the classical Buchstab function.
    """

    def _compute_steps(self, smax):
        self._steps = list(range(floor(smax) + 1))

    def _successeur(self, current_left):
        if current_left not in self._steps or current_left == self._steps[-1]:
            return None
        return current_left + 1

    def _centre(self, current_left):
        return current_left + 1 / 2

    def _predecessor_leq(self, vals, target):
        cand = [v for v in vals if v <= target]
        return max(cand) if cand else None

    def _adjust_constant(self, current_left):
        if current_left <= 1 or current_left == len(self._steps) - 1:
            return

        previous_left = current_left - 1

        val_left = self._approximants[previous_left](1)
        val_right = self._approximants[current_left](-1)

        correction = val_left - val_right
        self._approximants[current_left] += correction

    def _developpe(self, numerator, denominator, denominator_minorant):
        """
        Compute a Taylor approximation at 0, on the interval [-1, 1],
        of numerator / denominator, where numerator, denominator are
        polynomials with denominator linear (deg denominator = 1). It
        is assumed that denominator_minorant is a lower-bound for
        |denominator|. The method is borrowed from
        polynomial_ring. Bradshaw's sagemath implementation of the
        Dickman rho function.
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
        approximant = self._approximants[current_left]
        total_error = approximant.base_ring()(0)
        constant_error = approximant.coefficient(0).rad()
        truncation_degree = approximant.degree()
        while truncation_degree >= 0 and approximant.coefficient(truncation_degree).abs().upper() < constant_error:
            total_error += approximant.coefficient(truncation_degree).above_abs()
            truncation_degree -= 1
        self._approximants[current_left] = approximant.truncate(truncation_degree + 1) + total_error

    def _compute_approximants(self):
        for current_left in self._steps:
            if current_left in (0, 1, max(self._steps)):
                continue

            previous_left = current_left - 1

            previous_center = previous_left + 1 / 2

            polynomial_ring = self._polynomial_ring
            scalar_ring = self._scalar_ring
            poly_gen = polynomial_ring.gen()
            if previous_left in self._approximants:
                term = self._developpe(self._approximants[previous_left](poly_gen),
                                       poly_gen / 2 + scalar_ring(previous_center),
                                       scalar_ring(-1 / 2 + previous_center))
                self._approximants[current_left] = scalar_ring(1 / 2) * term.integral()
                self._adjust_constant(current_left)
                self._truncate(current_left)

    def __init__(self, precision=53):
        self._taylor_degree = floor(precision * RR(2).log() / RR(3).log())

        scalar_ring = RealBallField(precision)
        polynomial_ring = PolynomialRing(scalar_ring, "t")
        self._scalar_ring = scalar_ring
        self._polynomial_ring = polynomial_ring

        self._compute_steps(3)

        approximants = {}
        approximants[0] = polynomial_ring(0)
        approximants[1] = polynomial_ring(1)
        self._approximants = approximants

        self._compute_approximants()

    def __call__(self, sval):
        scalar_ring = self._scalar_ring
        if sval >= self._steps[-1]:
            self._compute_steps(RR(sval) * 1.1 + 1)
            self._compute_approximants()

        steps = self._steps

        if sval <= 0:
            return scalar_ring(0)

        current_left = self._predecessor_leq(steps, sval)
        if current_left is None:
            return scalar_ring(0)

        current_center = self._centre(current_left)

        return self._approximants[current_left](2 * (scalar_ring(sval) - scalar_ring(current_center)))

    def base_ring(self):
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
