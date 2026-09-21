# TODO

## Regarding ``number_theoretic_dde_solutions.py``

* Include the DickmanRho function to the existing class in ``sage.functions.transcendental`` (https://github.com/sagemath/sage/blob/develop/src/sage/functions/transcendental.py)
* Include the Buchstab omega or B function as well.
* Not sure about the Friedlander σ function.

## Regarding ``sieve_integral.py``

* The most proper way to integrate this to Sage would be to extend the ``integrate`` method from ``Polyhedron`` to handle piece-wise rational fractions with coefficients in a RealBallField.
* Let ``sieve_integral`` take an arbitrary integrand which admits polynomial approximations on a partition of the polytope: a factor would be described by the partition it is constant on, the approximant attached to each piece, and a bound on the approximation error, and the integral would be computed piece by piece as it already is for the Buchstab factor. Together with it, a way to build such a factor is needed, general enough to cover both what the library handles today, `B(t_r/t_s)` for the Harman examples, and the two-factor integrand of `h(α, β)`, `ρ(affine) ω(affine)`, whose factors are composed with affine forms rather than with a coordinate ratio. This would let ``Computation of h/computation_of_h.sage`` drop ``erreur_absolue`` and ``calcule_I_polyedre_equilibre``, which exist only because the library cannot express that integrand; note that ``h_bar2`` also does its own cross-term error bookkeeping around the call, which would have to move into the library.
