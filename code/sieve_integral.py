from functools import partial
from multiprocessing import Pool as pool
from operator import le, ge, lt, gt, eq
from sage.rings.real_arb import RealBallField
from sage.rings.rational_field import QQ
from sage.rings.real_mpfr import RealField, RR
from sage.functions.other import floor, ceil, binomial
from sage.numerical.mip import (MixedIntegerLinearProgram,
                                MIPSolverException)
from sage.matrix.special import diagonal_matrix
from sage.rings.power_series_ring import PowerSeriesRing
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.big_oh import O
from sage.misc.misc_c import prod
from sage.geometry.polyhedron.constructor import Polyhedron
from sage.geometry.polyhedron.base import Polyhedron_base
from sage.symbolic.expression import Expression
from tqdm import tqdm

load("number_theoretic_dde_solutions.py")

# Set the following to 1 if you do not wish to parallelize the
# computation.

NUM_PROCESSORS = 7


def expression_to_polynomial(expression, variables):
    r"""
    Convert a symbolic expression to a rational polynomial.

    INPUT:
    
    - ``expression`` -- symbolic expression

    - ``variables`` -- a list or tuple of variables containing all
    variables in `expression`

    OUTPUT: a tuple of

    - a rational polynomial,

    - generators of its polynomial ring, in the same order as the
    variables.

    EXAMPLES:

    sage: var("a, b, c")
    (a, b, c)

    sage: expression_to_polynomial(c^2 - a + 2*b + 4, (a, b, c))
    (x2^2 - x0 + 2*x1 + 4, (x0, x1, x2))

    sage: expression_to_polynomial(c^2 - a + 2*b + 4, (c, b, a))
    (x0^2 + 2*x1 - x2 + 4, (x0, x1, x2))

    sage: expression_to_polynomial(sin(a) + b, (a, b))
    TypeError     Traceback (most recent call last)
    ...
    TypeError: The expression to be converted does not coerce to a
    rational polynomial.

    sage: expression_to_polynomial(a^2 + b - c, (a, b))
    ValueError    Traceback (most recent call last)
    ...
    ValueError: Expression contains variables not in the variable list
    """


    if not set(variables).issuperset(set(expression.variables())):
        raise ValueError("Expression contains variables not in the "
                         "variable list")
    ring = PolynomialRing(QQ, len(variables), "x")
    generators = ring.gens()
    try:
        polynomial = (
            ring(expression.subs({variables[i]: generators[i]
                                  for i in range(len(variables))}))
        )
    except TypeError as exc:
        raise TypeError("The expression to be converted does not "
                        "coerce to a rational polynomial.")
    return (polynomial, generators)



def polynomial_to_equation(polynomial, generators):
    r"""Return the coefficients which represent the inequality
    ``polynomial >= 0``.

    INPUT:

    - ``polynomial``: a polynomial of degree at most 1.

    - ``generators``: a list of tuple of generators of the parent ring
    of ``polynomial``.

    OUTPUT: a tuple ``(a_0, a_1, ... a_k)`` where ``a_0`` is the
    constant coefficient of p, and ``a_i`` are the coefficient of the
    ``i``-th generator in ``polynomial``.


    EXAMPLES:

    sage: R.<x, y, z, t> = PolynomialRing(QQ)

    sage: generators = R.gens(); generators
    (x, y, z, t)

    sage: polynomial_to_equation(x + 3*y - z + 3, generators)
    (3, 1, 3, -1, 0)

    sage: polynomial_to_equation(x + y + z + t - 1, generators)
    (-1, 1, 1, 1, 1)

    sage: polynomial_to_equation(x^2 + y, generators)
    ValueError                   Traceback (most recent call last)
    ...
    ValueError: The polynomial should be affine linear.

    """

    if polynomial.degree() > 1:
        raise ValueError("The polynomial should be affine linear.")
    return ((polynomial.constant_coefficient(),)
            + tuple(polynomial.monomial_coefficient(x) for x in generators))





def symbolic_to_eqns(expressions):
    r"""
    Convert an iterable of symbolic linear equalities or equalities
    into a list of coefficients suitable for use in Polyhedron.

    INPUT:
    
    - ``expressions``: an iterable (list, tuple or set for instance)
    of symbolic linear equalities or inequalities.

    OUTPUT: a tuple ``(eqns, ieqs, variables)`` consisting of

    - ``eqns``: a list of equalities, where an equality
    `Ax + b = 0` is represented by a list ``(b, A)``

    - ``ieqs``: similarly, a list of inequalities, with
    `Ax + b \geq 0` represented by a list ``(b, A)``

    - ``variables``: a list of symbolic variables, ordered in such a
    way that the ``i``-th element of ``variables`` corresponds to the
    ``i``-th coordinate in the equalities and inequalities.

    EXAMPLES:

    sage: from operator import le, ge, lt, gt, eq

    sage: var("a, b, c")

    sage: symbolic_to_eqns({2*a + b - c <= 0, a > c, b + a == 3*c})
    ([(0, 1, -3, 1)], [(0, -2, 1, -1), (0, 1, -1, 0)], (a, c, b))

    sage: symbolic_to_eqns({a == 2*b, b == 2*a, a + b >= 0})
    ([(0, 1, -2), (0, -2, 1)], [(0, 1, 1)], (a, b))

    sage: symbolic_to_eqns({sqrt(2) * a >= b})
    TypeError                    Traceback (most recent call last)
    ...
    TypeError: The expression to be converted does not coerce
    to a rational polynomial.

    sage: symbolic_to_eqns({a^2 <= b, a == 0, b == 0})
    ValueError                   Traceback (most recent call last)
    ...
    ValueError: It is likely that one of the inequalities provided
    involves non-linear terms.
    
    .. SEEALSO::

    :meth:`sage.geometry.polyhedron.constructor.Polyhedron`

    """

    ieqs = []
    eqns = []
    variables = set()
    for expr in expressions:
        for x in expr.variables():
            variables.add(x)
    variables = list(variables)
    for expression in expressions:
        if not isinstance(expression, Expression):
            raise ValueError
        if expression.operator() not in (le, ge, lt, gt, eq):
            raise ValueError
        try:
            polynomial, generators = expression_to_polynomial(
                expression.lhs() - expression.rhs(),
                variables
            )
        except ValueError as exc:
            raise ValueError(
                "It is likely that one of the inequalities provided "
                "involves non-polynomial terms."
            ) from exc

        if expression.operator() in (le, lt):
            polynomial = -polynomial
        try:
            if expression.operator() == eq:
                eqns.append(polynomial_to_equation(polynomial, generators))
            else:
                ieqs.append(polynomial_to_equation(polynomial, generators))
        except ValueError as exc:
            raise ValueError(
                "It is likely that one of the inequalities provided "
                "involves non-linear terms."
            ) from exc
    return (eqns, ieqs, variables)



def are_inequalities_compatible(ieqs, border = False):
    r"""
    Checks if a linear system of linear inequalities has a solution.

    INPUT:

    - ``ieqs``: a list of inequalities `Ax + b \geq 0` in the shape
    ``(b, A)``.

    - ``border`` (default: ``False``): if ``True``, checks for the
    existence of a solution under the additional condition
    `\sum_i x_i = 1`.

    OUTPUT: True or False.

    TESTS:

    sage: from sage.numerical.mip import MIPSolverException

    sage: # Testing {0 ≤ x ≤ y and x + y ≤ 1}
    sage: ieqs = [(1, -1, -1), (0, -1, 1), (0, 1, 0)]
    sage: are_inequalities_compatible(ieqs)
    True

    sage: # Testing {2 ≤ x and 2 ≤ y and x + y ≤ 3}
    sage: ieqs = [(3, -1, -1), (-2, 1, 0), (-2, 0, 1)]
    sage: are_inequalities_compatible(ieqs)
    False

    sage: # Testing {x, y, z ≥ 1/2 and x + y + z = 1}
    sage: ieqs = [(-1/2, 1, 0, 0), (-1/2, 0, 1, 0), (-1/2, 0, 0, 1)]
    sage: are_inequalities_compatible(ieqs, border = True)
    False

    """

    if not ieqs:
        return True

    dim = len(ieqs[0]) - 1
    lin_prog = MixedIntegerLinearProgram()
    x = lin_prog.new_variable(real = True)
    if border:
        lin_prog.add_constraint(
            sum(x[j] for j in range(1, dim+1)) == 1
        )
    for ieq in ieqs:
        if len(ieq) != dim + 1:
            raise ValueError("The length of the inequalities "
                             "provided are inconsistent")
        lin_prog.add_constraint(
            ieq[0] + sum(x[j] * ieq[j] for j in range(1, dim+1)) >= 0
        )
        lin_prog.set_objective(None)
    try:
        lin_prog.solve(objective_only = True)
        return True
    except MIPSolverException:
        return False


def polytope_if_dim(ieqs, dim, border):
    r"""
    Check whether inequalities define a polytope of a given dimension.
    
    The expected dimension is ``dim`` if ``border = False``,
    ``dim - 1`` otherwise.

    INPUT:

    - ``ieqs``: a list of inequalities `Ax + b \geq 0` in the shape
    ``(b, A)``.

    - ``dim``: an integer.

    - ``border``: a boolean.

    OUTPUT: Either a polytope, or None, depending on whether the
    given inequalities produce a polyhedron of the dimension expected
    from the parameters ``dim`` and ``border``.


    EXAMPLES:

    sage: ieqs = [(1, -1, -1), (0, 1, 0), (0, -1, 1)]
    sage: polytope_if_dim(ieqs, 2, False)
    A 2-dimensional polyhedron in QQ^2 defined as the convex hull of
    3 vertices

    sage: polytope_if_dim(ieqs, 3, False)
    None

    sage: ieqs = [(1, -1, -1), (0, 1, -1), (0, -1, 1)]
    sage: polytope_if_dim(ieqs, 2, False)
    None
    """

    if are_inequalities_compatible(ieqs, border):
        polytope = Polyhedron(
            eqns = [(-1,) + (1,)*dim] if border else [],
            ieqs = ieqs
        )
        if polytope.dimension() == (dim - 1 if border else dim):
            return polytope
    return None


def absolute_error(polynomial):
    r"""
    Compute a midpoint polynomial approximation together with a
    uniform absolute error bound.

    The input polynomial is assumed to be univariate with
    coefficients in a ``RealBallField``, and its variable is assumed
    to vary in [-1, 1].

    INPUT:

    - ``polynomial`` -- a univariate polynomial over a
    ``RealBallField``.

    OUTPUT: a tuple ``(polynomial_mid, error)``, where
    
    - ``polynomial_mid`` is a polynomial over a ``RealField``
    (same precision) whose coefficients are the midpoints of those
    of ``polynomial``,
    - ``error`` is the sum of the radii of the coefficients
    of ``polynomial``.

    EXAMPLES:

    sage: K = RealBallField(30)
    sage: KP.<t> = PolynomialRing(R)
    sage: polynomial = sum(K(RR.random_element()) * t^j
                           for j in range(5))
    sage: x0 = random() * 2 - 1
    sage: 



    """
    scalar_field = polynomial.base_ring()
    error_field = RealField(prec=scalar_field.precision(), rnd="RNDU")
    poly_ring = PolynomialRing(QQ, "t")
    poly_gen = poly_ring.gen()
    coeffs = polynomial.monomial_coefficients()
    polynomial = poly_ring(sum(coeff.mid() * (poly_gen**degree)
                               for degree, coeff in coeffs.items()))
    erreur_abs = sum(error_field(coeff.rad()) for coeff in coeffs.values())
    return (polynomial, erreur_abs)


def make_ieq_fulldim_from_eqn(ieq, eqn):
    
    r"""
    Eliminate the last coordinate from an inequality
    using an equation.

    INPUT:

    - ``ieq`` -- a list representing an inequality
    `A_1 x + b_1 \geq 0`, in Polyhedron format ``(b1, A1)``.

    - ``eqn`` -- a list representing an equality
    `A_2 x + b_2 = 0`, also in Polyhedron format ``(b2, A2)``.
    The last entry must be nonzero

    OUTPUT:

    A list of length ``len(ieq) - 1`` representing the inequality
    obtained after eliminating the last coordinate.

    EXAMPLES::

        sage: make_ieq_fulldim_from_eqn((1, 2, 3), (4, 5, 6))
        (-1, -1/2)

    The last coefficient of ``eqn`` must be nonzero::

        sage: make_ieq_fulldim_from_eqn((1, 2, 3), (4, 5, 0))
        ValueError              Traceback (most recent call last):
        ...
        ValueError: Cannot eliminate the target variable. Check its
        coefficient in the equality

    TESTS::

        sage: make_ieq_fulldim_from_eqn((0, 0, 1), (1, -1, 2))
        (-1/2, 1/2)

    SEEALSO::

        :class:`~sage.geometry.polyhedron.constructor.Polyhedron`,
        for details on the format of equalities and inequalities.
    """

    dim = len(ieq)-1
    if eqn[dim].is_zero():
        raise ValueError("Cannot eliminate the target variable. "
                         "Check its coefficient in the equality")
    return tuple(ieq[i] - ieq[dim] * eqn[i] / eqn[dim]
                 for i in range(dim))


def latte_integrate(polytope, polynomial = None):
    """
    Call the LattE software to integrate polynomial on polytope.

    There should be at most one equation in the polytope's H representation.
    If this is the case, integration is performed with respect to the Lebesgue
    measure on all coordinates except the last.

    INPUT:

    - ``polytope`` -- the polytope to integrate on, as a Polyhedron
    object.

    - ``polynomial`` (default: ``None``) -- a multivariate rational polynomial.

    OUTPUT: the integral of polynomial over polytope, computed by LattE using
    the cone decomposition method.


    EXAMPLES:

    sage: polytope = Polyhedron([(0, 0), (0, 1), (1, 0)])
    sage: R.<t1, t2> = PolynomialRing(QQ)
    sage: latte_integrate(polytope, t1 * t2^2)
    1/60

    sage: polytope = Polyhedron(eqns = [(-3, 1, 2)],
                                ieqs = [(0, 1, 0), (0, 0, 1)])
    sage: latte_integrate(polytope, t1 * t2^2)
    27/16
    
    sage: polytope = polytopes.simplex(2)
    sage: R.<t1, t2, t3> = PolynomialRing(QQ)
    sage: latte_integrate(polytope, t1 * t2^2 * t3^3)
    1/3360

    sage: polytope = Polyhedron(eqns = [(-1, 1, 0)],
                                ieqs = [(0, 1, 0), (0, -1, 1), (2, 0, -1)])
    sage: R.<t1, t2> = PolynomialRing(QQ, 2)
    sage: latte_integrate(polytope, t1 * t2)
    0
    sage: polytope.integrate(t1 * t2, measure = 'induced')
    3/2

    sage: polytope = (polytopes.hypercube(5)
                      .faces(face_dimension=3)[0]
                      .as_polyhedron(base_ring = QQ))
    sage: polytope
    A 3-dimensional polyhedron in QQ^5 defined as the convex hull of 8 vertices
    
    sage: R.<t1, t2, t3, t4, t5> = PolynomialRing(QQ)
    sage: latte_integrate(polytope, (t1*t2*t3*t4*t5)^2)
    NotImplementedError                       Traceback (most recent call last)
    ...
    NotImplementedError: The latte_integrate helper function is not designed
    for polytopes with two defining equations or more. Use
    polytope.integrate(...) instead.
    
    sage: polytope.integrate((t1*t2*t3*t4*t5)^2, measure="induced")
    8/27

    """
    
    if len(polytope.equations_list()) not in (0, 1):
        # We only handle the case of a singe equality
        # of type sum_j r_j x_j = 1
        raise NotImplementedError(
            "The latte_integrate helper function is not designed for "
            "polytopes with two defining equations or more. "
            "Use polytope.integrate(...) instead."
        )

    if not polytope.equations_list():
        # We are integrating over a full-dimensional polytope
        if polynomial is None:
            return polytope.volume()
        return polytope.integrate(
            polynomial.change_ring(QQ),
            algorithm = 'cone-decompose'
        )

    original_ieqs = polytope.inequalities_list()

    eqn = polytope.equations_list()[0]
    if eqn[-1].is_zero():
        return 0

    dim = polytope.ambient_dimension()
    if polynomial is None:
        old_ring = PolynomialRing(QQ, "x", dim)
        poly = old_ring.one()
    else:
        poly = polynomial.change_ring(QQ)
        old_ring = poly.parent()
    old_gens = old_ring.gens()
    new_ring = PolynomialRing(QQ, "t", dim - 1)
    new_gens = new_ring.gens()
    dic = {}
    for i in range(dim - 1):
        dic[old_gens[i]] = new_gens[i]
    dic[old_gens[dim-1]] = (
        -(eqn[0] + sum(new_gens[i] * eqn[1+i]
                       for i in range(dim - 1))) / eqn[dim]
    )
    polynomial_new = new_ring(poly.subs(dic))

    polytope_new = Polyhedron(
        ieqs = [make_ieq_fulldim_from_eqn(ieq, eqn)
                for ieq in original_ieqs]
    )

    return polytope_new.integrate(polynomial_new,
                                  algorithm = 'cone-decompose')



class PolytopeSummary(SageObject):
    r"""
        Carries a list of inequalities defining a polytope, along with some
        geometric data.

    It is assumed that the defining polytope's equations list is either
        empty, or else containes the single equation `\sum_i x_i = 1`.

    .. WARNING::

    This class was designed for internal use in sieve_integral only, and
        used in situations where the polytope is either full-dimensional,
        or else cut out by the specific equation `\sum_i x_i = 1`.
        The constructor checks whether the number of equations is zero or
        one, and thereafter discards the input polytope's equations list.
        It *DOES NOT* check that the defining equation is the one above, and
        will therefore yield wrong results if it is called with a polytope
        cut out by one equation, other than `\sum_i x_i = 1`.


    
    EXAMPLES:

    sage: polytope = polytopes.hypercube(3); polytope
        A 3-dimensional polyhedron in ZZ^3 defined as the convex hull
        of 8 vertices
        sage: summary = PolytopeSummary(polytope)
        sage: summary.unpack()
        ([[1, -1, 0, 0],
        ...
        [1, 0, 1, 0]],
        False,
        (-1, -1, -1),
        (1, 1, 1),
        8)

    sage: polytope = polytopes.simplex(2); polytope
        A 2-dimensional polyhedron in ZZ^3 defined as the convex hull
        of 3 vertices
        sage: summary = PolytopeSummary(polytope)
        sage: summary.unpack()
        ([[1, 0, -1, -1], [0, 0, 1, 0], [0, 0, 0, 1]],
        True,
        (0, 0, 0),
        (1, 1, 1),
        1/2)

    sage: polytope = (polytopes.hypercube(5)
                      .faces(face_dimension=3)[0]
                      .as_polyhedron(base_ring = QQ))
        sage: polytope
        A 3-dimensional polyhedron in QQ^5 defined as the convex hull of
        8 vertices
        sage: summary = PolytopeSummary(polytope)
        NotImplementedError               Traceback (most recent call last)
        ...
        NotImplementedError: The polytope should be full-dimensional, or
    else cut out by the single equation (-1, 1, ..., 1).



    """

    def __init__(self, polytope):
        r"""
        INPUT:

        - ``polytope`` -- a polytope as a Polyhedron object.

        """

        if len(polytope.equations_list()) not in (0, 1):
            raise NotImplementedError(
                "The polytope should be full-dimensional, or "
                "else cut out by the single equation (-1, 1, ..., 1)."
            )
        self._border = (len(polytope.equations_list()) == 1)
        self._ieqs = polytope.inequalities_list()
        self._mini, self._maxi = polytope.bounding_box()
        self._volume = polytope.volume(engine = 'latte',
                                       algorithm = 'cone-decompose',
                                       measure = 'induced_rational')

    def _repr_(self):
        return (
            "Polytope summary data "
            + f"with {len(self._ieqs)} inequalities, "
            + ("cut out by the equation (-1, 1, ..., 1), " if self._border
               else "full-dimensional, ")
            + f"and volume {self._volume}"
        )

    def volume(self):
        return self._volume

    def unpack(self):
        r"""
        Return all the data contained in this polytope summary.

        OUTPUT:

        - ``ieqs`` -- the original polytope's inequalities.

        - ``border`` -- whether the polytope is cut out
        by the equation `\sum_i x_i = 1`.

        - ``mini`` -- the lower coordinates of a bounding box, as a tuple.
        
        - ``maxi`` -- the upper coordinates of a bounding box, as a tuple.

        - ``volume`` -- the polytope's volume, either as a full-dimensional
        polytope (if ``border`` is ``False``) or else projected onto any
        ``d-1`` of the coordinates (if ``border`` is ``True``, with ``d``
        the dimension of the polytope).

        """

        return (self._ieqs, 
                self._border,
                self._mini,
                self._maxi,
                self._volume)



def eqn_sum_one(dim):
    r"""
    Return the inequality `\sum_{i=1}^d x_i = 1` in the Polyhedron
    ``(b, A)`` format.

    INPUT:

    - ``dim`` -- the dimension

    OUTPUT:

    - the tuple ``(-1, 1, ..., 1)``, where ``1`` is repeated ``dim`` times.

    EXAMPLES:

    sage: eqn_sum_one(4)
    (-1, 1, 1, 1, 1)

    sage: polytope = Polyhedron(eqns = [eqn_sum_one(2)])
    sage: polytope
    A 1-dimensional polyhedron in QQ^2 defined as the convex hull of
    1 vertex and 1 line

    sage: polytope.equations()
    (An equation (1, 1) x - 1 == 0,)

    sage: Polyhedron(eqns = [eqn_sum_one(0)])
    The empty polyhedron in QQ^0

    sage: eqn_sum_one(-2)
    ValueError                        Traceback (most recent call last)
    ...
    ValueError: The dimension should be non-negative.

    """

    if dim<0:
        raise ValueError("The dimension should be non-negative.")
    return (-1,) + (1,) * dim


def balanced_polytope_integrate(polytope_summary,
                                truncation_degree,
                                scalar_field,
                                verbose = 0):
    r"""
    Approximate the integral over a polytope of the Taylor expansion
    main term of `1/x_1 ... x_d`.

    The integrand `1/x_1 ... x_d` is Taylor-expanded at the midpoint of the
    polytope, and truncated to degree ``truncation_degree - 1``.
    
    .. WARNING::
    The Taylor truncation error is not estimated in this method. This method
    approximates only the integral of the main term in the Taylor expansion.

    INPUT:

    - ``polytope_summary`` (a PolytopeSummary instance)
       -- a data containing the polytope inequalities, the coordinates of
    a bounding box, and its volume.

    - ``truncation_degree`` -- the degree at which we truncate the Taylor
    expansion of ``1/x_1 ... x_k`` at the mid-points of the polytope.

    - ``scalar_field`` -- the RealBallField in which the value will be returned.

    - verbose (integer) -- if positive, will show debugging information.

    OUTPUT:
    an element of ``scalar_field`` containing the value
    `\int_P dt / x_1 ... x_k` if `P` is full-dimensional, where `P` is
    the polytope. If `P` should be cut out by the equation `\sum_i x_i = 1`,
    then the integral is computed with respect to `dx_1 ... dx_{d-1}` with `d`
    the ambient dimension.

    EXAMPLES:

    sage: polytope = polytopes.hypercube(3) + vector((5, 5, 5))
    sage: sage: polytope.bounding_box()
    ((4, 4, 4), (6, 6, 6))

    sage: balanced_polytope_integrate(PolytopeSummary(polytope),
                                      truncation_degree = 4,
                                      scalar_field = RBF).endpoints()
    (0.0665599999999999, 0.0665600000000001)

    sage: approximant = taylor(1/(x*y*z), (x, 5), (y, 5), (z, 5), 2)
    sage: polynomial = approximant.polynomial(QQ)
    sage: polytope.integrate(polynomial)
    208/3125

    sage: RR(208/3125)
    0.0665600000000000

    sage: RR(log(6/4)^3)
    0.0666592560084858

    sage: balanced_polytope_integrate(PolytopeSummary(polytope),
    ....:                             truncation_degree = 23,
    ....:                             scalar_field = RBF).endpoints()
    (0.0666592560084857, 0.0666592560084858)


    sage: polytope = Polyhedron(eqns = [eqn_sum_one(2)],
                                ieqs = [(-2/5, 1, 0), (-2/5, 0, 1)])
    sage: polytope.bounding_box()
    ((2/5, 2/5), (3/5, 3/5))

    sage: summary = PolytopeSummary(polytope)
    sage: balanced_polytope_integrate(summary,
    ....:                             truncation_degree = 3,
    ....:                             scalar_field = RBF).endpoints()
    (0.810666666666666, 0.810666666666667)

    sage: approximant = taylor(1/(x*y), (x, 1/2), (y, 1/2), 2)
    sage: polynomial = approximant.polynomial(QQ)
    sage: RR(polytope.integrate(polynomial, measure = "induced") / sqrt(2))
    0.810666666666667

    sage: exact_value = integrate(1 / (x*(1-x)), (x, 2/5, 3/5))
    sage: exact_value
    2*log(3/5) - 2*log(2/5)
    sage: RR(exact_value)
    0.810930216216329

    sage: balanced_polytope_integrate(summary,
    ....:                             truncation_degree = 25,
    ....:                             scalar_field = RBF).endpoints()
    (0.810930216216328, 0.810930216216329)

    """

    if truncation_degree <= 0:
        return scalar_field.zero()
    printifdbg = print if verbose > 0 else (lambda *x:None)
    precision = scalar_field.precision()

    ieqs, border, mini, maxi, vol = polytope_summary.unpack()
    if len(ieqs) == 0:
        return scalar_field.zero()
    polytope = Polyhedron(
        ieqs = ieqs,
        eqns = [eqn_sum_one(len(ieqs[0])-1)] if border else [],
        base_ring = QQ
    )

    majorant_product = 1/prod(mini)
    trivial_bound = vol * majorant_product

    if trivial_bound < 2**(-precision):
        trivial_bound_lower = vol / prod(maxi)
        valeur = (trivial_bound + trivial_bound_lower) / 2
        estimation_queue = (trivial_bound - trivial_bound_lower) / 2
        printifdbg(f"  I ≃ {RR(valeur)}, "
                   f"tail ≤ {RR(estimation_queue)} (trivial)")
        return scalar_field(valeur).add_error(estimation_queue)

    dim = polytope.ambient_dimension()
    if dim <= 1:
        raise NotImplementedError(
            "The polytope should be of dimension 2 or more. "
            "If you are integrating over an interval, use symbolic "
            "computation instead."
        )
    centers = tuple((mini[j]+maxi[j])/2 for j in range(dim))
    rayons = tuple((maxi[j]-mini[j])/2 for j in range(dim))
    ratios = tuple(rayons[j] / centers[j] for j in range(dim))
    translated_polytope = (
        polytope
        .translation(-c for c in centers)
        .linear_transformation(
            diagonal_matrix(1/ray for ray in rayons)
        )
    )
    dilation_factor = prod(rayons[:-1]) if border else prod(rayons)
    # We make it so that the new variables are between -1 and 1

    series_gens = PowerSeriesRing(QQ, "x", num_gens = dim).gens()

    fonction = (prod((1 + series_gens[j] * ratios[j]
                      + O(series_gens[j] ** truncation_degree))
                     .inverse() for j in range(dim))
                .polynomial())
    homog = fonction.homogeneous_components()

    truncation_degree = truncation_degree-1
    estimation_queue = 0
    for degree in range(truncation_degree, -1, -1):
        if degree in homog:
            upper_bound = sum(
                abs(coeff) for coeff in (
                    homog[degree].monomial_coefficients().values()
                )
            )
            new_majo = (
                scalar_field(vol * upper_bound / prod(centers))
                .upper()
            )
            if estimation_queue + new_majo > 2**(-precision):
                break
            estimation_queue += new_majo
            truncation_degree = degree

    polynome_latte = sum(
        homog_part
        for d, homog_part in homog.items()
        if d <= truncation_degree
    )
    valeur = (
        latte_integrate(translated_polytope,
                        polynome_latte)
        * dilation_factor / prod(centers)
    )
    printifdbg(f"  I ≃ {RR(valeur)}, tail ≤ {RR(estimation_queue)}"
               f" (degree {truncation_degree})")
    return scalar_field(valeur).add_error(estimation_queue)



def errorbound_inverse_tail(coeff, rho, dim, scalar_field):
    r"""
    Compute an upper-bound on the Taylor truncation of the series
    `1/(1 - x_1) \dotsb (1 - x_d)` to a prescribed precision.

    INPUT:

    - ``coeff``, ``rho`` -- rational numbers.

    - ``dim`` -- an integer.

    - ``scalar_field`` -- a RealBallField.

    OUTPUT: a tuple consisting of two objects:

    - ``degree``, a non-negative integer,

    - ``epsilon``, a RealField element with upper-rounding error,

    in such a way that :
    * the Taylor truncation error for the rational fraction
    `1/(1 - x_1) \dotsb (1 - x_d)` on `[-\rho, \rho]^{dim}`
    at degree ``degree`` is at most ``epsilon`` in absolute values,
    * the value ``coeff * coeff`` is bounded by the absolute precision
    of ``scalar_field``.

    
    EXAMPLES:

    sage: errorbound_inverse_tail(1, 1/3, 1, RBF)
    (34, 8.99432546225716e-17)

    sage: RR((1/3)^34 / (1 - 1/3))
    8.99432546225715e-17

    sage: degree, epsilon = errorbound_inverse_tail(7, 1/4, 3, RBF)
    sage: (degree, epsilon)
    (33, 1.09644964147099e-17)

    sage: S.<u, v, w> = PowerSeriesRing(QQ)
    sage: series = ((1+u)*(1+v)*(1+w) + O(u, v, w)^33).inverse()
    sage: approximant = series.polynomial()
    sage: exact_error = (4/3)^3 - approximant(-1/4, -1/4, -1/4)
             # The error is attained at u, v, w minimal
    sage: RR(exact_error)
    1.09644964147099e-17

    sage: 7 * epsilon < 2^(-RBF.precision())
    True

    """
    if dim <= 0:
        raise ValueError("The dimension should be a positive integer.")
    if rho <= 0:
        raise ValueError("The Taylor expansion range is empty.")
    if rho >= 1:
        raise ValueError("The Taylor expansion range contains zero.")

    precision = scalar_field.precision()
    ratio = rho/(1-rho)
    def relative_error(degree):
        # Regarding the following bound, see Lemma 2.2 of "Computing
        # sieve integrals using LattE, and the density of integers
        # with a localized divisor" (Drappeau, Mounier)
        return (
            (rho ** (degree-1))
            * sum((ratio**(dim-j)) * binomial(degree+dim-1, j)
                  for j in range(dim))
        )
    degree = 0
    while coeff * relative_error(degree) > 2**(-precision):
        degree += 1
    epsilon = scalar_field(relative_error(degree)).upper()
    return (degree, epsilon)


def is_polytope_on_sum1(polytope):
    r"""
    Check if the polytope is of the hyperplane `\sum_i x_i = 1`, or if it
    is full-dimensional.

    Raises an exception if none of the two cases occurs.

    INPUT: ``polytope`` (Polyhedron object) -- a polytope

    OUTPUT: a boolean.

    Returns True if ``polytope`` is in the hyperplane `\sum_i x_i = 1` and
    is of full-dimensional inside it, and False if ``polytope`` is full-
    dimensional.


    EXAMPLES:

    sage: is_polytope_on_sum1(polytopes.hypercube(4))
    False

    sage: is_polytope_on_sum1(polytopes.simplex(4))
    True

    sage: is_polytope_on_sum1(polytope_line)
    ValueError                                Traceback (most recent call last)
    ...
    ValueError: This method is currently handling only polytopes which are
    either full-dimensional, or cut out by the specific equation
    `\sum_i x_i = 1`.

    """

    if not isinstance(polytope, Polyhedron_base):
        raise TypeError("The polytope should be a Polyhedron object.")

    if len(polytope.equations_list()) not in (0, 1):
        raise ValueError(
            "This method is currently handling only polytopes which are "
            "either full-dimensional, or cut out by the specific equation "
            r"`\sum_i x_i = 1`."
        )

    if len(polytope.equations_list()) == 0:
        return False

    eqn = polytope.equations_list()[0]
    if all(x == -eqn[0] for x in eqn[1:]):
        return True
    raise ValueError(
        "This method is currently handling only polytopes which are "
        "either full-dimensional, or cut out by the specific equation "
        r"`\sum_i x_i = 1`."
    )


def balance_polytope(polytope, facteur, verbose = 0):
    r"""
    Partition a polytope into pieces where the coordinates vary within a
    certain factor of each other.

    INPUT:

    - ``polytope`` (Polyhedron object) -- a polytope.

    - ``facteur`` -- a number greater than 1.

    - ``verbose`` -- a boolean.

    OUTPUT: a list of elements of class PolytopeSummary, each of which
    describes a polytope `Q` which is balanced in the sense that for every
    `x, y \in Q`, we have `x_i / y_i \leq facteur` for every `i`.


    EXAMPLES:

    sage: polytope = polytopes.hypercube(3) + vector((4, 4, 4))
    sage: balanced_polytopes = balance_polytope(polytope, 1.3)
    sage: balanced_polytopes
    [Polytope summary data with ... volume 729/1000,
     Polytope summary data with ... volume 891/1000,
     Polytope summary data with ... volume 891/1000,
     Polytope summary data with ... volume 1089/1000,
     Polytope summary data with ... volume 891/1000,
     Polytope summary data with ... volume 1089/1000,
     Polytope summary data with ... volume 1089/1000,
     Polytope summary data with ... volume 1331/1000]

    sage: add(Q.volume() for Q in balanced_polytopes)
    8

    sage: polytope = vector((1/10,) * 4) + 6/10 * polytopes.simplex(3)
    sage: polytope.vertices()
    (A vertex at (1/10, 1/10, 1/10, 7/10),
     A vertex at (1/10, 1/10, 7/10, 1/10),
     A vertex at (1/10, 7/10, 1/10, 1/10),
     A vertex at (7/10, 1/10, 1/10, 1/10))
    sage: balanced_polytopes = balance_polytope(polytope, 1.5)

    sage: balanced_polytopes
    [Polytope summary data with 6 inequalities, cut out by the equation
      (-1, 1, ..., 1), and volume 1/1000,
    ...
     Polytope summary data with 4 inequalities, cut out by the equation
      (-1, 1, ..., 1), and volume 1/6000]

    sage: add(Q.volume() for Q in balanced_polytopes)
    9/250

    sage: polytope.volume(measure='induced_rational')
    9/250

    """
    printifdbg = print if verbose > 0 else (lambda *x:None)
    dim = polytope.ambient_dimension()

    border = is_polytope_on_sum1(polytope)
    
    facteur = RR(facteur)
    if facteur <= 1:
        raise ValueError("The balancing parameter must be greater than 1.")

    def facteur_pow(i):
        return QQ(round(facteur**i, 3))
    # TODO à affiner

    # First compute the number of cuts to be made in each dimensions.
    nombre_decoupes = []
    mini, maxi = polytope.bounding_box()
    nombre_decoupes = tuple(
        ceil((maxi[j]/mini[j]).log() / facteur.log()) for j in range(dim)
    )

    printifdbg(f"Number of slices = {nombre_decoupes}")

    # Compute the inequalities lists which describe the cuts, in each
    # dimension.
    planches = []
    for j in range(dim):
        ieq_coeff_plus = (0,) * j + (1,) + (0,) * (dim-1-j)
        ieq_coeff_minus = (0,) * j + (-1,) + (0,) * (dim-1-j)
        planches.append(tuple(
            [(-mini[j]*facteur_pow(exponent),) + ieq_coeff_plus,
             (mini[j]*facteur_pow(exponent+1),) + ieq_coeff_minus]
            for exponent in range(nombre_decoupes[j])
        ))

    ieqs_pile = [polytope.inequalities_list()]

    printifdbg("Slicing...")
    for dim_decoupe in range(dim):
        ieqs_pile_new = []
        for ieqs_old in ieqs_pile:
            for decoupe in planches[dim_decoupe]:
                ieqs_candidate = ieqs_old + decoupe
                if are_inequalities_compatible(ieqs_candidate):
                    ieqs_pile_new.append(ieqs_candidate)
        ieqs_pile = list(ieqs_pile_new)

    printifdbg("Computing volumes and bounding box...")

    split_polytopes = []
    for ieqs in ieqs_pile if verbose <= 0 else tqdm(ieqs_pile):
        current_polytope = polytope_if_dim(ieqs, dim, border)
        if current_polytope is not None:
            split_polytopes.append(
                PolytopeSummary(current_polytope)
            )

    printifdbg(f"Slicing done, there are {len(split_polytopes)} polytopes.")

    return split_polytopes


def process_polytope(polytope_summary,
                     scalar_field,
                     verbose = 0):
    """This is a helper function which processes the sieve integral
    over a balanced element Qel of the partition of our original
    polytope. This function is for use in multiprocessing.

    """
    border, mini, maxi, vol = polytope_summary.unpack()[1:]
    dim = len(mini)
    # Compute the largest ratio between values, in each dimension,
    # and use this value to estimate how far we need to truncate
    facteur = max(maxi[j]/mini[j] for j in range(dim))
    truncation_degree, prod_rel_error = errorbound_inverse_tail(
        vol/prod(mini), (facteur-1)/(facteur+1), dim, scalar_field
    )
    val = balanced_polytope_integrate(polytope_summary,
                                      truncation_degree,
                                      scalar_field,
                                      verbose = verbose)
    majorant_product = 1 / prod(mini)
    abs_error = vol * majorant_product * prod_rel_error
    return scalar_field(val).add_error(abs_error)



def sieve_integral(polytope_data, precision = 20, facteur = 1.3, verbose = 0):
    r"""
    Compute the integral over the given polytope of
    ``dt_1 ... dt_d / t_1 ... t_d``.

    This is currently implemented only in the case if the polytope described
    by ``polytope_data`` is full-dimensional, or if it is cut out by the
    equation `\sum_i x_i = 1`.


    INPUT:

    - ``polytope_data`` -- data describing the polytope we integrate over,
    either a genuine Polyhedron object, or an iterable of symbolic linear
    inequalities.

    - ``precision`` (default: 20 bits) -- the precision, in bits, at which the
    computation is carried out. Note that the final precision may deteriorate
    by a quantity depending on the number of balanced polytopes
    in the decomposition.

    - ``facteur`` -- a number greater than 1.

    - ``verbose`` (default: False) -- if positive, will print updates on the
    computation.


    OUTPUT: a ``RealBall`` containing the value of the requested integral.

    EXAMPLES:

    sage: polytope = Polyhedron([(1,), (2,)])
    sage: sieve_integral(polytope)
    [0.69315 +/- 4.37e-6]

    sage: polytope = (vector((2, 2, 2))
    ....:             + polytopes.hypercube(3, intervals='zero_one'))
    sage: polytope.bounding_box()
    ((2, 2, 2), (3, 3, 3))

    sage: sieve_integral(polytope)
    [0.06666 +/- 6.17e-6]

    sage: RR(log(3/2)^3)
    0.0666592560084858

    sage: polytope = Polyhedron(eqns = [(-1, 1, 1, 1)],
    ....:                       ieqs = [(-1/7, 1, 0, 0),
    ....:                               (0, -1, 1, 0),
    ....:                               (0, 0, -1, 1)])
    sage: sieve_integral(polytope)
    [0.9569 +/- 3.06e-5]

    sage: sieve_integral(polytope, precision = 50)
    [0.9569135271017 +/- 5.30e-14]

    sage: RR(1/2 * integrate(log(6 - x)/x, (x, 1, 5)))
    0.956913527101737

    """

    if not isinstance(polytope_data, Polyhedron_base):
        eqns, ieqs = symbolic_to_eqns(polytope_data)[:2]
        polytope_data = Polyhedron(ieqs = ieqs, eqns = eqns)

    polytope = Polyhedron(polytope_data, base_ring = QQ)

    printifdbg = print if verbose > 0 else (lambda *x:None)
    scalar_field = RealBallField(precision)
    if polytope.is_empty():
        return scalar_field.zero()

    border = is_polytope_on_sum1(polytope_data)

    mini, maxi = polytope.bounding_box()
    border = is_polytope_on_sum1(polytope)
    if polytope.ambient_dimension() == 1:
        if border:
            return scalar_field.one()
        return scalar_field(maxi[0] / mini[0]).log()
    
    majorant_product = 1/prod(mini)
    vol = polytope.volume(engine = 'latte',
                          algorithm = 'cone-decompose',
                          measure = 'induced_rational')
    trivial_bound = vol * majorant_product
    if trivial_bound < 2**(-precision):
        printifdbg(f"sieve_integral returning a trivial bound "
                   f"≤ {trivial_bound}")
        return scalar_field(trivial_bound/2, rad = trivial_bound/2)

    balanced_polytopes = balance_polytope(polytope,
                                          facteur,
                                          verbose = verbose - 1)
    process_polytope_partial = partial(process_polytope,
                                       scalar_field = scalar_field,
                                       verbose = verbose - 1)

    # Parallelize only if there are sufficiently many polytopes; there
    # is a bit of overhead in setting up the multiprocessing.
    if len(balanced_polytopes)>20:
        with pool(NUM_PROCESSORS) as p:
            reslist = p.map(process_polytope_partial,
                            balanced_polytopes)
    else:
        reslist = [process_polytope_partial(balanced_polytope)
                   for balanced_polytope in balanced_polytopes]
    sum_values = sum(sorted(reslist))
    printifdbg(f"sieve_integral returning {sum_values}")
    return sum_values



def polynomial_trivial_bound(polynomial, maxi):
    r"""
    Compute a trivial bound on ``polynomial`` under the hypothesis that the
    variables are bounded above in absolute values by ``maxi``.

    INPUT:

    - ``polynomial`` -- a polynomial with coefficients in a RealBallField.

    - ``maxi`` -- a list or tuple of numbers coercible to RealBallField.

    OUTPUT:

    - a RealField element with upper-rounding.

    EXAMPLE:

    sage: K = PolynomialRing(RBF, "x")
    sage: polynomial = K(chebyshev_T(3, x))
    sage: polynomial
    4.000000000000000*x^3 - 3.000000000000000*x

    sage: polynomial_trivial_bound(polynomial, (1,))
    7.00000000000000

    sage: expansion = taylor(1/(1 + x + y + z),
    ....:                    (x, 0), (y, 0), (z, 0), 10)
    sage: polynomial = expansion.polynomial(RBF)
    sage: polynomial_trivial_bound(polynomial, (1/4, 1/4, 1/4))
    3.83105945587159

    """

    coeffs = polynomial.monomial_coefficients()
    if polynomial.parent().ngens() == 1:
        return sum(coeff.above_abs().upper() * maxi[0] ** d
                   for d, coeff in coeffs.items())
    return sum(coeff.above_abs().upper() * prod(maxi[j] ** d[j]
                                                for j in range(len(d)))
               for d, coeff in coeffs.items())




def sieve_integral_harman(polytope_data,
                          buchstab_variables,
                          integrand_polynomial_data = None,
                          precision = 20,
                          verbose = 0):
    r"""
    Compute the integral over a polytope of the function given by
    `B(x_i/x_j) F(x) / x_1 \dotsb x_d`, where `B(u) = u \omega(u)`
    with `\omega` the Buchstab function, and `F` is a polynomial.

    INPUT:

    - ``polytope_data`` -- a Polyhedron object, or an iterable of symbolic
    equality and/or inequalities.

    - ``buchstab_variables`` -- a tuple of variables, or indices, describing
    the argument of the Buchstab function in the integrand.

    - ``integrand_polynomial_data`` -- either a symbolic expression in
    variables which appear in ``polytope_data``, or a Polynomial.

    - ``precision`` (default: 20) -- the requested precision, in bits.

    - ``verbose`` (default: 0) -- if positive, gives information on the
    computation.

    OUTPUT: a RealBall containing the value of the integral.

    EXAMPLES:
    sage: polytope = {x + y == 1, 1/4 < x, x < y}
    sage: polynomial = 1 + x
    sage: sieve_integral_harman(polytope, (y, x), polynomial)
    [1.6921 +/- 2.88e-5]

    sage: A = integrate((1+x) / (x*(1-x)), (x, 1/3, 1/2))
    sage: B = integrate((1+log(1/x-2)) * (1+x) / (x*(1-x)), (x, 1/4, 1/3))
    sage: RR(A + B)
    1.69211856327702

    sage: sieve_integral_harman(polytope, (y, x), polynomial, precision=50)
    [1.6921185632770 +/- 3.27e-14]

    """

    if not isinstance(polytope_data, Polyhedron_base):
        eqns, ieqs, variables = symbolic_to_eqns(polytope_data)
        polytope_data = Polyhedron(ieqs = ieqs, eqns = eqns)
    else:
        variables = None

    polytope = Polyhedron(polytope_data, base_ring = QQ)
    if polytope.is_empty():
        return RealBallField(precision=precision).zero()
    border = is_polytope_on_sum1(polytope)
    dim = polytope.ambient_dimension()

    if buchstab_variables is None:
        buchstab_indices = None
    elif len(buchstab_variables) != 2:
        raise ValueError(
            "The buchstab_variables parameter should be a tuple or "
            "list with two elements, either two symbolic variables, "
            "or two integers."
        )
    else:
        bv_numerator, bv_denominator = buchstab_variables
        if variables is not None:
            if (bv_numerator not in variables
                or bv_denominator not in variables):
                raise ValueError(
                    "Each of the Buchstab variables should be one "
                    "of the symbolic variables appearing "
                    "in the inequalities"
                )
            buchstab_indices = (variables.index(bv_numerator),
                                variables.index(bv_denominator))
        else:
            buchstab_indices = (int(bv_numerator), int(bv_numerator))

    if integrand_polynomial_data is None:
        integrand_polynomial = None
    else:
        poly_ring = PolynomialRing(QQ, "x", len(variables))
        poly_gens = poly_ring.gens()
        if isinstance(integrand_polynomial_data, Expression):
            try:
                integrand_polynomial = poly_ring(
                    integrand_polynomial_data.subs(
                        {variables[j]: poly_gens[j]
                         for j in range(len(variables))}
                    )
                )
            except TypeError as exc:
                raise (
                    TypeError("The integrand symbolic expression "
                              "does not parse to a polynomial.")
                ) from exc
        else:
            try:
                integrand_polynomial = (
                    poly_ring(integrand_polynomial_data)
                )
            except TypeError as exc:
                raise TypeError("The integrand polynomial cannot "
                                "be parsed. Check that its base ring "
                                "coerces to QQ, and that it has the "
                                "right number of variables.") from exc

    # We ignore all possible equalities implemented in P, and
    # implement the border only.  Beware that if your polytope had to
    # complementary inequalities, Sage might merge them into an
    # equality, which would disappear here.
    if polytope.is_empty():
        return 0
    ieqs = polytope.inequalities_list()
    dim = polytope.ambient_dimension()
    if border:
        eqns = [eqn_sum_one(dim)]
    else:
        eqns = None


    # We load an approximation of the Buchstab B(u) = u ω(u) function
    # from a specific module.
    buchstab_function = BuchstabB(precision = precision)
    scalar_field = buchstab_function.base_ring()
    mini, maxi = polytope.bounding_box()
    buchstab_function(
        maxi[buchstab_indices[0]] / mini[buchstab_indices[1]]
    )  # TODO à affiner

    buchstab_approximants_data = (
        buchstab_function.approximants_list()
    )
    # Now buchstab_approximants_data is a list of triplets (n, n+1,
    # P_n) where n is an integer, and P_n is a polynomial with
    # coefficients in a RealBallField, such that B(u) = P_n(2*u-2*n-1)
    # for every u in (n, n+1). Note that the argument of P_n is always
    # a number between -1 and 1.

    value = scalar_field.zero()

    for (
        left_endpoint,
        right_endpoint,
        buchstab_approximant
    ) in buchstab_approximants_data:
        ieq_left = [0] * (dim+1)
        ieq_left[buchstab_indices[0]+1] = 1
        ieq_left[buchstab_indices[1]+1] = -left_endpoint
        ieq_right = [0] * (dim+1)
        ieq_right[buchstab_indices[0]+1] = -1
        ieq_right[buchstab_indices[1]+1] = right_endpoint
        # The inequalities ieq_left and ieq_left implement the
        # restriction "t_r / t_s in (left_endpoint, right_endpoint)".
        polytope = Polyhedron(
            eqns = eqns,
            ieqs = ieqs + [tuple(ieq_left), tuple(ieq_right)]
        )
        if (
            not polytope.is_empty()
            and (not border and polytope.dimension() == dim
                 or border and polytope.dimension() == dim-1)
        ):
            value += sieve_integral_polyratio(
                polytope,
                buchstab_approximant,
                QQ((left_endpoint + right_endpoint)/2),
                buchstab_indices,
                polynomial2 = integrand_polynomial,
                verbose = verbose,
                precision = precision
            )


    return value


def sieve_integral_polyratio(polytope,
                             polynomial1,
                             midpt,
                             ratio_idx,
                             polynomial2 = None,
                             verbose = 0,
                             precision = 20):
    """
    Computes the sieve integral over the Polyhedron P, of
    polynomial1(2*(t_r / t_s - midpt)) * polynomial2(t) / prod(t_j).
    Here (r, s) is ratio_idx.  It is assumed that P is such that the
    argument of polynomial1 here is of absolute value at most 1.

    A more detailed docstring will be provided once we figure out how best
    to integrate this into the sage libraries.
    """
    printifdbg = print if verbose > 0 else (lambda *x:None)
    scalar_field = RealBallField(precision)

    border = is_polytope_on_sum1(polytope)
    mini, maxi = polytope.bounding_box()
    majorant_polynomial1 = (
        1 if polynomial1 is None
        else polynomial_trivial_bound(polynomial1, (1,))
    )
    majorant_polynomial2 = (
        1 if polynomial2 is None
        else polynomial_trivial_bound(polynomial2.change_ring(scalar_field),
                                      maxi)
    )
    majorant_integrant = majorant_polynomial1 * majorant_polynomial2
    majorant_product = 1/prod(mini)
    trivial_bound = (polytope.volume(measure = 'induced_rational')
                     * majorant_product * majorant_integrant)

    if trivial_bound < 2**(-precision):
        printifdbg(f"sieve_integral_polyratio returning the "
                   f"trivial bound ≤ {trivial_bound}")
        return scalar_field(trivial_bound/2, rad = trivial_bound/2)

    split_polytopes = balance_polytope(polytope,
                                       facteur = 1.3,
                                       verbose = verbose - 1)
    process_polytope_partial = partial(process_polytope_polyratio,
                                       scalar_field = scalar_field,
                                       polynomial1 = polynomial1,
                                       midpt = midpt,
                                       polynomial2 = polynomial2,
                                       ratio_idx = ratio_idx,
                                       verbose = verbose - 1)
    if len(split_polytopes) > 20:
        with pool(NUM_PROCESSORS) as p:
            values = p.map(process_polytope_partial, split_polytopes)
    else:
        values = [process_polytope_partial(polytope_summary)
                  for polytope_summary in split_polytopes]
    sum_values = sum(sorted(values))
    printifdbg(f"sieve_integral_polyratio returning {sum_values}")
    return sum_values

def process_polytope_polyratio(polytope_summary,
                               scalar_field,
                               polynomial1,
                               midpt,
                               ratio_idx,
                               polynomial2 = None,
                               verbose = 0):
    """
    This is a helper function which mainly calls the balanced
    integration routine on the polytope described by Qel.  Qel is a
    PolytopeSummary Computes the sieve integral over the polytope
    P, of polynomial1(2*(t_r / t_s - midpt)) * polynomial2(t) /
    prod(t_j).  Here (r, s) is ratio_idx.  It is assumed that P is
    such that the argument of polynomial1 here is of absolute value at
    most 1.

    A more detailed docstring will be provided once we figure out how best
    to integrate this into the sage libraries.
    """

    border, mini, maxi, vol = polytope_summary.unpack()[1:]
    dim = len(mini)
    facteur = max(maxi[j]/mini[j] for j in range(dim))
    truncation_degree, prod_rel_error = (
        errorbound_inverse_tail(vol/prod(mini),
                                (facteur-1)/(facteur+1),
                                dim,
                                scalar_field)
    )
    value = balanced_polytope_integrate_polyratio(
        polytope_summary,
        polynomial1,
        midpt,
        ratio_idx,
        polynomial2,
        truncation_degree,
        scalar_field,
        verbose = verbose - 1
    )
    majorant_product = 1/prod(mini)
    majorant_polynomial1 = (
        1 if polynomial1 is None
        else polynomial_trivial_bound(polynomial1, (1,))
    )
    majorant_polynomial2 = (
        1 if polynomial2 is None
        else polynomial_trivial_bound(polynomial2.change_ring(scalar_field),
                                      maxi)
    )
    majorant_integrant = majorant_polynomial1 * majorant_polynomial2
    abs_error = (
        vol * majorant_integrant * majorant_product * prod_rel_error
    )
    return scalar_field(value).add_error(abs_error)



def balanced_polytope_integrate_polyratio(
        polytope_summary,
        polynomial1,
        midpt,
        ratio_idx,
        polynomial2,
        truncation_degree,
        scalar_field,
        verbose = 0
):
    """
    Computes the sieve integral over the polytope P, of
    polynomial1(2*(t_r / t_s - midpt)) * polynomial2(t) * Π(t), where
    Π(t) is the order N truncation of 1/prod(t_j), and (r, s) is
    ratio_idx.  It is assumed that P is such that the argument of
    polynomial1 here is of absolute value at most 1.

    A more detailed docstring will be provided once we figure out how best
    to integrate this into the sage libraries.
    """
    if truncation_degree == 0:
        return scalar_field.zero()
    printifdbg = print if verbose > 0 else (lambda *x: None)
    precision = scalar_field.precision()

    if polynomial1 is not None:
        poly1, poly1_epsilon = absolute_error(polynomial1)
    else:
        poly1_epsilon = scalar_field.zero().upper()
        poly1 = None

    ieqs, border, mini, maxi, vol = polytope_summary.unpack()

    majorant_polynomial1 = (
        1 if poly1 is None
        else polynomial_trivial_bound(polynomial1, (1,))
    )
    majorant_polynomial2 = (
        1 if polynomial2 is None
        else polynomial_trivial_bound(polynomial2.change_ring(scalar_field),
                                      maxi)
    )
    majorant_integrant = majorant_polynomial1 * majorant_polynomial2
    majorant_product = 1/prod(mini)
    trivial_bound = vol * majorant_product * majorant_integrant

    if trivial_bound < 2**(-precision):
        printifdbg(f"  |I| ≤ {RR(trivial_bound)}")
        # Without additional information the polynomials could be
        # negative.
        return scalar_field.zero().add_error(trivial_bound)

    if len(ieqs) == 0:
        return scalar_field.zero()
    dim = len(ieqs[0]) - 1
    if dim <= 1:
        raise NotImplementedError(
            "The polytope should be of dimension 2 or more. "
            "If you are integrating over an interval, use symbolic "
            "computation instead."
        )
    polytope = Polyhedron(ieqs = ieqs,
                          eqns = [eqn_sum_one(dim)] if border else [],
                          base_ring = QQ)

    centers = tuple((mini[j] + maxi[j]) / 2 for j in range(dim))
    rayons = tuple((maxi[j] - mini[j]) / 2 for j in range(dim))
    ratios = tuple(rayons[j] / centers[j] for j in range(dim))
    translated_polytope = (
        polytope
        .translation(-c for c in centers)
        .linear_transformation(
            diagonal_matrix(1 / ray for ray in rayons)
        )
    )
    dilation_factor = prod(rayons[:-1]) if border else prod(rayons)
    # In case border = True, we will project onto the first dim-1
    # coordinates before integrating.

    series_ring = PowerSeriesRing(QQ, "x", num_gens = dim)
    series_gens = series_ring.gens()
    fonction = prod(
        (1 + series_gens[j] * ratios[j]
         + O(series_gens[j]**truncation_degree)).inverse()
        for j in range(dim)
    ).polynomial()

    if poly1 is not None:
        numer_idx, denom_idx = ratio_idx
        truncation_degree_denom = (
            1 + floor(precision * RR(2).log()
                      / (1/ratios[denom_idx]).log())
        )
        denom = (1 + series_gens[denom_idx] * ratios[denom_idx]
                 + O(series_gens[denom_idx]**truncation_degree_denom))
        inv_denom = denom.inverse()
        numer = (
            centers[numer_idx]
            + series_gens[numer_idx] * rayons[numer_idx]
        )
        poly1_arg = (
            (numer * inv_denom / centers[denom_idx] - midpt) * 2
        )
        poly1_arg_max = maxi[numer_idx] / mini[denom_idx]
        poly1_arg_eps = (
            ratios[denom_idx]**truncation_degree_denom * poly1_arg_max
        )
        poly1_mc = poly1.monomial_coefficients()
        poly1_val_eps = (
            poly1_epsilon
            + sum(coeff * ((poly1_arg_max + poly1_arg_eps)**d
                           - poly1_arg_max**d)
                  for d, coeff in poly1_mc.items())
        )
        poly1_compose = poly1.subs(poly1_arg).polynomial()
    else:
        poly1_compose = 1
        poly1_val_eps = scalar_field.zero()

    if polynomial2 is None:
        poly2_compose = series_ring.one().polynomial()
    else:
        p2vars = polynomial2.parent().gens()
        poly2_compose = series_ring(
            polynomial2.subs({
                p2vars[j]: centers[j] + series_gens[j] * rayons[j]
                for j in range(dim)
            })
        ).polynomial()
    homog = (
        fonction
        * (poly1_compose if poly1 is not None else 1)
        * poly2_compose
    ).homogeneous_components()


    dmax = truncation_degree - 1
    estimation_queue = scalar_field(
        vol * majorant_product * majorant_polynomial2 * poly1_val_eps
    ).upper()
    for degree in range(dmax, -1, -1):
        if degree in homog:
            upper_bound = (
                sum(abs(coeff) for coeff in (homog[degree]
                                             .monomial_coefficients()
                                             .values()))
            )
            new_majo = (
                scalar_field(vol * upper_bound / prod(centers))
                .upper()
            )
            if estimation_queue + new_majo > 2**(-precision):
                break
            estimation_queue += new_majo
            dmax = degree


    polynomial_latte = sum(homog_part
                           for degree, homog_part
                           in homog.items()
                           if degree <= dmax)
    valeur = (
        latte_integrate(translated_polytope, polynomial_latte)
        * dilation_factor / prod(centers)
    )
    printifdbg(f"  I ≃ {RR(valeur)}, tail ≤ {RR(estimation_queue)} "
               f"(degree {dmax})")

    return scalar_field(valeur).add_error(estimation_queue)
