# Computation of h(α, β)

This folder contains the notebook used for the computation of $`h(\alpha, \beta)`$
described in *Computing sieve integrals using LattE, and the density of integers with a localized divisor* (Drappeau, Mounier).

* The Sage notebook ``Computation of h(alpha, beta).ipynb`` contains the code relative to the computation.
* The script ``computation_of_h.sage`` is the same code, cell for cell, runnable outside Jupyter:
``cd "Computation of h" && sage computation_of_h.sage`` computes the three worked values of the notebook,
and ``load("computation_of_h.sage")`` in a Sage session started here defines everything and computes nothing.
Its docstring lists the places where it departs from the notebook.
* The HTML file ``Plot of h(alpha, beta).html`` contains the 3d surface plot containing our estimate graph of $`h(\alpha, \beta)`$
* The Dickman, Buchstab and Friedlander classes are loaded from ``../code/number_theoretic_dde_solutions.py``
* The script ``IncreasingBooleanFunctions.py`` contains the computations relative to Boolean functions which is used in the notebook.
This will be replaced by the use of Boolean antichains, once we figure out how to pass to the dual.
* The folder ``data`` contains raw data used for the Boolean functions script.
* The NumPy archive ``list_infos_ab.npz`` contains, for each k, the inequalities of the projections of the
polytopes Π_F to the (α, β) plane, which ``whichF`` uses in the notebook to find the F for which Δ_F(α, β) is non-empty.
It is written by ``ecrit_liste_ab`` and read by ``lit_liste_ab``.
* The JSON file ``hbar-echantillon.json`` contains the sample of values of $`\bar{h}(\alpha, \beta)`$ used for the plot,
as returned by ``h_bar``: α and β as exact rationals, and the real balls as their midpoint and radius, both exact
rationals, together with the precision. It is written by ``ecrit_liste`` and read by ``lit_liste``.

## Where changes go

From now on, changes to this computation go into ``computation_of_h.sage``. The notebook is kept as it is,
as the record of the published outputs: it is the only committed trace of the values of $`h(\alpha, \beta)`$
and of the figures made from them.

