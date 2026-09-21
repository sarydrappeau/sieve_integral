# Fixes

A record of the defects found in this code and of what was done about them, one
row per finding, with the evidence that established it. The delay-differential
equation module is covered here; later work extends the file.

Everything below was checked against SageMath 10.7 and LattE 1.7.6, and
`cd code && sage -t --long number_theoretic_dde_solutions.py` passes.

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

## Effects on the recorded values

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
