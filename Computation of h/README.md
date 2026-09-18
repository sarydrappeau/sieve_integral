# Computation of h(α, β)

This folder contains the notebook used for the computation of $`h(\alpha, \beta)`$
described in *Computing sieve integrals using LattE, and the density of integers with a localized divisor* (Drappeau, Mounier).

* The Sage notebook ``Computation of h(alpha, beta).ipynb`` contains the code relative to the computation.
* The HTML file ``Plot of h(alpha, beta).html`` contains the 3d surface plot containing our estimate graph of $`h(\alpha, \beta)`$
* The script ``IncreasingBooleanFunctions.py`` contains the computations relative to Boolean functions which is used in the notebook.
This will be replaced by the use of Boolean antichains, once we figure out how to pass to the dual.
* The folder ``data`` contains raw data used for the Boolean functions script.
* The NumPy archive ``list_infos_ab.npz`` contains, for each k, the inequalities of the projections of the
polytopes Π_F to the (α, β) plane, which ``whichF`` uses in the notebook to find the F for which Δ_F(α, β) is non-empty.
It is written by ``ecrit_liste_ab`` and read by ``lit_liste_ab``.

