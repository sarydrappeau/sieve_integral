# TODO

## Regarding ``number_theoretic_dde_solutions.py``

* Include the DickmanRho function to the existing class in ``sage.functions.transcendental`` (https://github.com/sagemath/sage/blob/develop/src/sage/functions/transcendental.py)
* Include the Buchstab omega or B function as well.
* Not sure about the Friedlander σ function.

## Regarding ``sieve_integral.py``

* The most proper way to integrate this to Sage would be to extend the ``integrate`` method from ``Polyhedron`` to handle piece-wise rational fractions with coefficients in a RealBallField.
