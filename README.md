# Computing of sieve integrals using LattE

Python / Sage script files providing methods to compute integrals of the shape

$$ I(P) = \int_P \frac{dt_1 \dotsb dt_k}{t_1 \dotsb t_k} $$

where $`P \subset {\mathbb R}_{>0}^k`$ is a rational convex polytope, and variants such as

$$ I^*(P) = \int_{P \cap H_d} \frac{dt_1 \dotsb dt_{k-1}}{t_1 \dotsb t_k} $$

where $`H = \lbrace (t_1, \dotsc, t_k) \in {\mathbb R}^k, \sum_i t_i = 1\rbrace`$.

## Description

### ``sieve_integral.py``

The script ``sieve_integral.py`` defines the methods ``sieve_integral(...)`` and ``sieve_integral_harman(...)``, among others.

* The method ``sieve_integral`` computes the integral of the shape $`I(P)`$ or $`I^*(P)`$.

* The method ``sieve_integral_harman`` computes integrals of the shape

$$ \int_P B\Big(\frac{t_i}{t_j}\Big) f(t_1, \dotsc, t_k) \frac{ dt_1 \dotsb dt_k}{t_1 \dotsb t_k}, $$

  where $`B(u) = u \omega(u)`$ and $`f`$ is a polynomial, as well as their variant

$$ \int_{P \cap H_k} B\Big(\frac{t_i}{t_j}\Big) f(t_1, \dotsc, t_k) \frac{ dt_1 \dotsb dt_{k-1}}{t_1 \dotsb t_k} $$

  where $`H_k = \lbrace(t_1, \dotsc, t_k)\in {\mathbb R}^k, \sum_i t_i = 1\rbrace`$.


### ``number_theoretic_dde_solutions.py``

The script ``number_theoretic_dde_solutions.py`` defines three classes:

- ``DickmanRho``: defines the [Dickman ρ function](https://en.wikipedia.org/wiki/Dickman_function).
- ``BuchstabB``: defines the function $u \mapsto u \omega(u)$, where $\omega$ is the [Buchstab ω function](https://en.wikipedia.org/wiki/Buchstab_function).
- ``FriedlanderS``: defines a variant of the [Friedlander σ function](https://doi.org/10.1112/plms/s3-33.3.565).

Compared with the current SageMath implementation of the Dickman Rho function, their specificity is to carry out the computation in a RealBallField.
The algorithm is the same as in the SageMath ``DickmanRho`` implementation, based on the [Marsaglia-Zaman-Marsaglia method](https://doi.org/10.2307/2008355).

### Examples

The folder ``Examples`` contains the notebook describing the script's output on a handful of cases from the sieve literature.
See the README.md file for more information about those.


## Getting Started

### Dependencies

* Requires Sage and LattE (see [Getting Started](#getting-started)).
* It was designed with Sage 10.7 and LattE 1.7.6. It may well work with earlier versions. Please report any incompatibility issue you witness.

### Installing

- Install SageMath following the instructions described here : [Sage Installation Guide](https://doc.sagemath.org/html/en/installation/index.html).
- Install ``LattE`` : [latte-integrale from conda-forge](https://anaconda.org/channels/conda-forge/packages/latte-integrale/overview)

### Executing program

- Download and place both scripts ``number_theoretic_dde_solutions.py`` and ``sieve_integral.py`` your favourite directory.
- Import or load them in your notebook using your favourite method.

## Help

The doc-strings contain informations relative to usage.
Feel free to contact [@sarydrappeau-arch](https://github.com/sarydrappeau-arch) if needed.

## Authors

* [Sary Drappeau](https://github.com/sarydrappeau-arch)
* Adrien Mounier

## Version History

* 0.1 (2026-06) Initial upload

## License

This project is licensed under the [GNU GPL 3 License](LICENSE).

## Acknowledgments

This work is based on the preprint *Computing sieve integrals using LattE, and the density of integers with a localized divisor* (Drappeau, Mounier).

Supported by CNRS, Université d'Aix-Marseille, Université Clermont-Auvergne, Institut Universitaire de France, ANR-20-CE91-0006 / FWF-I-4945-N, ANR-24-CE93-0016 / FNS-10.003.145.
