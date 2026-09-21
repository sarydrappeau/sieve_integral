# Fixes

A record of the defects found in this code and of what was done about them, one
row per finding, with the evidence that established it. The two modules of
`code/` are covered here; later work extends the file.

Everything below was checked against SageMath 10.7 and LattE 1.7.6, and
`cd code && sage -t --long number_theoretic_dde_solutions.py sieve_integral.py`
passes, with 57 and 174 tests.

## `code/number_theoretic_dde_solutions.py`

| # | Issue | Evidence | Resolution |
| --- | --- | --- | --- |
| 1 | `BuchstabB._truncate` accumulated the discarded coefficients as a value, so the approximant was moved by their sum where it had to be widened by it, and the ball could miss the true value | at 30 bits the values at `u = 13/2` and `15/2` do not overlap those of a 200-bit instance; interval by interval, the shift starts on `[2, 3]` at every precision and is of the order of the radius there (`1.9e-6` against `4.3e-6` at precision 20), outgrowing it as the intervals accumulate | accumulate into the radius, through `add_error`, as `DickmanRho` did |
| 2 | `FriedlanderS._truncate` was called from inside `if current_left == 1`, so one interval only was ever truncated | at 30 bits with `α = 2/5` the degrees stayed at 18 or 19, where `DickmanRho`'s fall from 15 to 0; the enclosure held, so this cost width and time rather than rigour | the call de-indented, the jump at 1 left inside the test |
| 3 | `FriedlanderS` accepted a parameter it cannot handle | `FriedlanderS(3/2)` left the interval at 1 without an approximant, since it has no predecessor at distance `α`: `S(5/4)` raised `KeyError` while `S(7/4)` quietly returned `1.0`; a non-rational `α` failed deep inside the construction, on a polynomial ring and a real field | `ValueError` unless `0 < α < 1`, `TypeError` if `QQ` cannot convert |
| 4 | `FriedlanderS.__call__` evaluated the polynomial at its argument without coercing it | `FriedlanderS(2/5)(RR(0.75))` and the same on a Python float raised `TypeError: no common canonical parent`, where `DickmanRho` and `BuchstabB` accept both | coerce into the ball field first |
| 5 | the grid of steps was built in floating point | `DickmanRho._centre(3)` returned the Python float `3.5`, and `_steps` held Python `int`s; `n + 1/2` and `±1/2` are exact in binary, so no value was wrong | the grid built in `QQ`, so that steps, centres and radii are exact |
| 6 | three near-identical copies of the same machinery, which had drifted apart | the drift is what produced rows 1 and 2: the truncation was right in one copy, shifted in another, and called once in the third | one base class, `DifferentialDelaySolution`; a subclass describes its own equation |
| 7 | extending the grid recomputed every interval, not only the new ones | it does not widen anything: an instance extended to `u = 8`, `9`, `5` has every coefficient identical, midpoint and radius, to a fresh instance built to the same reach, and six successive extensions give exactly the same thing as one | `_steps_to_compute`, the steps with no approximant yet |
| 8 | `Computation of h` held a second copy of the module, to be kept in sync by hand | it had drifted: it predated row 1, so the h computation, and the three h anchors of the test suite, still ran on the shifted approximants | deleted; the notebook, the script and the READMEs load `../code/number_theoretic_dde_solutions.py` |
| 9 | no doctests at all | — | 57, including the enclosure check across precisions that catches row 1 for all three classes at once |
| 10 | `_truncate`'s docstring claimed the discarded terms were bounded by the error in the constant coefficient | the loop compares each of them to that error separately and does not cap their total, so the radius it adds can come out a few times larger; Sage's own `dickman_rho` drops the same terms (`truncate_abs`) and, being floating point, does not account for them at all | the docstring corrected; the criterion left as it is |

### Effects on the recorded values

Row 1 moves every value computed through the Buchstab chain, which is 72 of the
80 published examples and the h computation. Comparing the reference values
before and after, exactly: 63 of the 83 entries moved, and every new ball
overlaps the entry it replaces, so no published value had been outside its
enclosure. The 20 that did not move are the eight examples which carry no
Buchstab factor, the nine whose ratio stays at or below 2 -- below `u = 2` the
approximants are the exact seeds, so nothing there can move -- and the three h
anchors, which were still being computed from the copy of row 8. Radii mostly widen, by a factor up to about 10, but can
also shrink, a wider constant term raising the threshold `_truncate` compares
against; the two largest changes are `stadlmann/42`, `9.08e-6` to `4.15e-5`, and
`stadlmann/34`, `9.68e-5` to `1.24e-4`.

The three h anchors keep their midpoints and only the last bits of one radius
move. `h_bar1` takes the difference of two values of `B` inside a single
interval, on which row 1 amounted to a constant shift, and a constant cancels in
a difference.

The two values recorded in the docstring of `sieve_integral_harman` were
re-measured: `2.88e-5` to `2.90e-5` at precision 20, and `3.27e-14` to `3.03e-14`
at precision 50, the second one tighter for the reason just given. Both still
contain `1.692118563277022181104155555999314...`.

`FriedlanderS` is wider than it was, row 2 having put its truncation back into
service: at `u = 17/5` and 30 bits, `2.36e-7` becomes `7.90e-6`. Nothing in the
repository consumes it.

## `code/sieve_integral.py`

Step 7 of the reorganisation. The first six rows are the ones that could make
the returned ball miss the true value.

| # | Issue | Evidence | Resolution |
| --- | --- | --- | --- |
| 1 | `balance_polytope` took its number of cuts from `ceil(RR log(maxi/mini) / RR log facteur)`, the smallest `n` with `facteur**n ≥ maxi/mini`; the cells are built from `QQ(round(facteur**i, 3))`, which rounds to nearest, so the last line could fall short of the box and drop a sliver | on `[1, 285605/100000] × [1, 6/5]`, where `round(1.3^4, 3) = 2.856` sits below `1.3^4 = 2.8561`, the pieces added to `232/625` against a volume of `37121/100000`, and `sieve_integral` returned `[0.19133226 +/- 2.47e-9]`, which does not contain `log(285605/100000)·log(6/5) = 0.19133545…` | the number of cuts is the smallest `n ≥ 1` with `mini_j · facteur_pow(n) ≥ maxi_j`, compared exactly in `QQ` on the values the cells are built from; no logarithm is involved |
| 2 | the truncation loop of `balanced_polytope_integrate_polyratio` started at `truncation_degree - 1`, the degree of the truncated series alone, where what it truncates is that series times the Buchstab approximant and the integrand polynomial; everything above the starting degree was discarded without being compared against the budget | on `merikoski/F2` at precision 30 the product reaches degree 20 where the loop started at 9, and the terms between weigh `3.11e-9` against a budget of `9.31e-10` | the loop starts at `max(homog)`, as the h notebook always did |
| 3 | `poly1_val_eps` summed the coefficients of the approximant with their signs, where it bounds an error | `7.21e-6` instead of `8.49e-6` on the first piece of the `sieve_integral_harman` docstring example, a shortfall larger than the `9.54e-7` budget on its own | `abs(coeff)` |
| 4 | `poly1_arg_eps`, the error made by truncating `1/(1 + x_s r_s)` below degree `N`, was short by a factor 2: the tail `r_s^N/(1 - r_s)` is carried through `t_r/c_s ≤ maxi_r/c_s` **and** the outer factor 2 of the argument `2(t_r/t_s − midpt)` | — | `2·maxi_r·r_s^N/mini_s`, using `c_s(1 − r_s) = mini_s` |
| 5 | `absolute_error` promised the midpoints of the coefficients but built its polynomial over `QQ` by coercion, and coercing a `RealNumber` into `QQ` reconstructs the simplest nearby rational; the displacement, up to half an ulp per coefficient, was covered by nothing, the error returned beside it being the sum of the radii | at 30 bits the midpoint of `K(1/3)` is `357913941/1073741824`, and the coercion gives `1/3`, `1.59e-7` away where half an ulp is `2.38e-7` | `coeff.mid().exact_rational()`; LattE is also 47% faster on five examples, dyadics sharing one power of two as a denominator where reconstructed fractions have coprime ones |
| 6 | `poly1_arg_max` was the raw ratio `maxi_r/mini_s`, which is not a bound on the argument of the approximant, and was raised to its degree | up to 3.2 on the published examples, where the argument is bounded by 1 | the caller intersects the polytope with `n < t_r/t_s < n+1` and passes `midpt = n + 1/2`, so the argument lies in `[-1, 1]`: `poly1_arg_max = 1` |
| 7 | `sieve_integral_harman` with a `Polyhedron` and integer Buchstab indices wrote both slice inequalities to the same slot, `(int(bv_numerator), int(bv_numerator))` | every slice came out empty and the function returned zero, where the symbolic form gives `[1.2458 +/- 4.17e-5]` | `int(bv_denominator)` for the second |
| 8 | the same function with a `Polyhedron` and an integrand sized its polynomial ring with `len(variables)`, where `variables` is `None` | `TypeError: object of type 'NoneType' has no len()` | the ring is sized from the ambient dimension, and a symbolic integrand with a `Polyhedron` is rejected with a message |
| 9 | `buchstab_variables=None` was documented and subscripted `None` | `TypeError: 'NoneType' object is not subscriptable` | it integrates `F(t)/t_1⋯t_d`, which `sieve_integral` cannot express: on `[1,2]^2` with `F = 1 + x_1` the result contains `(1 + log 2) log 2` |
| 10 | nothing checked that the polytope lies in the open positive orthant, which is a hypothesis of the method: the balancing takes logarithms of the coordinate ratios, and the expansions are centred at positive midpoints | `ZeroDivisionError: rational division by zero`, or, on the unit cube, `unable to convert ceil(3.81149468670840*I*pi) to an integer` | `check_positive_orthant` at both entry points |
| 11 | `latte_integrate` returned 0 when the last coefficient of the equation vanishes, that is when the coordinate it eliminates is not one the equation constrains; the value is right for the measure the docstring promises, the projection having collapsed to a point, but it answers a degenerate question with a plausible number | the induced integral of `t1·t2` over that segment is `3/2` | `NotImplementedError`, naming the two ways out: permute the variables, which makes the elimination possible and gives `3/2`, or integrate with `measure='induced'` |
| 12 | `symbolic_to_eqns` took the coordinate order from a `set`, an implementation detail, and its two `ValueError`s had no message and could not be reached, the loop collecting the variables calling `expr.variables()` before the checks | the order was `[c, b, a]` and `[t3, t2, t1]`; a string in the argument raised `AttributeError: 'str' object has no attribute 'variables'` | variables sorted by name; the collecting and the validating are one loop, and the messages name the expression at fault |
| 13 | `sage -t code/sieve_integral.py` aborted while *parsing*, so not one example had ever run | `ValueError: line 56 of the docstring for sieve_integral.PolytopeSummary has inconsistent leading whitespace`; behind it, continuation lines without `....:`, exception examples written as IPython sessions, symbolic examples with no `var()`, an example stopping mid-sentence, and a name the docstring never defined | 174 doctests, with `# optional - latte_int` on the blocks that need LattE |
| 14 | the file was a script, not a module: `load("number_theoretic_dde_solutions.py")` at top level, `SageObject` never imported, and `round` taken from the session | in a bare Python process `round(RealNumber, 3)` raises `type RealNumber doesn't define __round__`, the Sage session's `round` being `sage.misc.functional.round` | a module docstring, the GPL header, explicit imports, and the companion imported when it is on `sys.path` and `load`ed otherwise |
| 15 | `facteur` was a parameter of `sieve_integral` but hard-coded to 1.3 inside `sieve_integral_polyratio`, so the Buchstab chain ignored it; `NUM_PROCESSORS` was a hard-coded 7; the objective of the feasibility LP was set once per inequality; both workers unpacked a flag they do not use; and `sieve_integral_polyratio` measured its polytope without naming the engine | — | `facteur` and `processes` threaded through, `NUM_PROCESSORS = max(1, (os.cpu_count() or 2) - 1)`, one `sum_over_pieces` helper for the two mappings, and the rest cleaned |

### Effect on the recorded values

Rows 2 to 6 move 63 of the 83 references. Every new ball contains the one it
replaces, and every one is wider: the radius grows by a factor 1.07 in the
median, by 1.36 at most (`merikoski/F4`, `3.70e-15` to `5.03e-15`), and by more
than five per cent in 50 of the 63. The 20 that do not move are the examples
with no Buchstab factor, the seven Stadlmann polytopes whose approximant is the
constant `B = 1` of `[1, 2]`, and the three h anchors, which go through the
notebook's own integrator.

Row 1 moves nothing: over the 80 published examples, `balance_polytope` is
called 285 times, on 1276 axes, and not one grid failed to reach its bounding
box; the new count agrees with the old on every one of them. Row 12 moves 11
Stadlmann examples, ten of them by one unit in the last place of the radius --
same pieces, same midpoint, a different summation order -- and `stadlmann/02`
from `[0.11600 +/- 7.71e-6]` to `[0.11328 +/- 6.92e-6]`, which the old reference
does not contain: that polytope comes out with two pieces instead of three,
because the floating-point LP rejects feasible cells and which ones depends on
the conditioning, hence on the coordinate order. Both values are too small; the
one with the cells decided exactly is about 0.1221.

The suite takes 2227s where it took 1761s, +26%, every run parallel and on an
idle machine. Most of it is row 2, the only change which must cost time, since
it keeps terms that used to be dropped: measured right after it the suite took
2035s, three fifths of the whole increase, and isolating the changes over five
examples puts all of that on the honest truncation start, row 5 being a saving.
The rest follows the coordinate order of row 12, which decides the coordinate
`latte_integrate` eliminates and so the polytope LattE decomposes: at ten bits
above its own precision `stadlmann/17` takes 34.1s sorted against 29.6s in the
order a set gave, while `maynard/I6` and `ford-maynard/I5` went the other way,
by 13% and 11%. The micro-benchmarks, which cover the plain chain and the
differential-delay solutions rather than the Buchstab chain, do not move: 7.53s
against 7.39s, a ratio of 0.98.

### Left for later

* `Computation of h/computation_of_h.sage` has the same defect as row 5 in its
  own chain, and this step does not repair it: `erreur_absolue` builds its
  midpoint polynomial over `Reals(prec)`, where the midpoints are exact, but the
  polynomial reaches LattE through `latte_integrate`, whose `change_ring(QQ)` is
  the same reconstruction, again covered by nothing. On `h/0.45-0.5` at
  precision `(20, 14)` all 43 polynomials handed to LattE have their
  coefficients in `RealField(14)`, the largest carrying 220 of them, and the
  reconstruction moves them by up to `3.0e-5` each and `2.8e-4` summed over one
  polynomial, against the `2^-14 = 6.1e-5` that computation asks for. Its integrator is not the
  library's, and cannot be until `sieve_integral` can express a product of
  interval-wise polynomials composed with affine forms -- the entry `TODO.md`
  now carries.
* `stadlmann/02` is still the floating-point LP above, which step 9 decides
  exactly.

## Traps worth remembering

* `QQ` on a real *reconstructs* a nearby rational by continued fractions rather
  than taking the exact binary value: `QQ(0.1)` is `1/10`, and `QQ(RR(pi))` is
  `245850922/78256779`. `RealNumber.exact_rational()` is the exact conversion.
* `RealBall.add_error` on a ball argument uses that ball's *upper* bound.
* The module is `load`ed as a plain `.py` file, which the Sage preparser does not
  see, so `1 / 2` in it is a Python float.
* Importing a Sage submodule, such as `sage.functions.other`, before `sage.all`
  in a fresh process trips a circular import inside Sage itself:
  `AttributeError: cannot access submodule 'function' of module 'sage.symbolic'`.
  Everything is taken from `sage.all`.
* Far from the first intervals the accumulated radius outgrows the value: at 30
  bits `ρ(12)` is only known to be below `1e-5`, and comparing `ρ(30) ≈ 3.27e-50`
  with Sage's `dickman_rho` needs 300 bits to say anything.
* In a Sage session `round` is `sage.misc.functional.round`, which shadows the
  builtin and accepts a `RealNumber`. Code meant to be imported has to say so.
  The two are not interchangeable either: `round(x, 3)` and
  `(x*1000).round()/1000` disagree on a tie, and give different grids at
  `facteur = 1.5`.
* A doctest whose polytope exceeds twenty pieces fails. The pieces then go to a
  worker pool, which pickles its worker by name, and under `sage -t` the
  module's functions are not reachable as `__main__.process_polytope`.
* Which cells the feasibility LP wrongly rejects depends on the conditioning,
  so a value it affects moves when the coordinates are reordered.
* `PowerSeriesRing(QQ, "x", num_gens=d)` truncates by *total* degree, which is
  why the plain chain has no counterpart to row 2: its expansion stops exactly
  at `truncation_degree - 1`.
