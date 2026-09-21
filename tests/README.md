# Tests

A small suite for the sieve integrals, run with Sage: the scripts it
exercises need the Sage preparser and `load`. **Run everything from
the repository root.**

```
sage tests/run_tests.sage                  # the quick tests, about 20 seconds
sage tests/run_tests.sage --long           # everything, expect about thirty minutes
sage tests/run_tests.sage --only-timing    # the micro-benchmarks
sage tests/run_tests.sage 'regression/chen*' 'closedform/*'   # globs on names
sage tests/run_tests.sage --long --list    # print the selected names, run nothing
```

The runner exits non-zero if a test fails. A test listed in `KNOWN_ISSUES` is
expected to fail: its failure is reported as `xfail` and does not fail the run,
and if it passes, the line reads `XPASS`, which is worth looking into.

## What is checked

`tests/test_examples.sage` covers the 80 examples of `Examples/` and the three
h tests of `Computation of h/computation_of_h.sage`:

| name | what it checks |
| --- | --- |
| `regression/<example>` | the value overlaps the one in `reference_values.json` |
| `selfconsistency/<example>` | the value overlaps the value at ten more bits |
| `bound/stadlmann/NN` | the value is below the bound of Stadlmann's paper |
| `notebook/h/...` | the h test overlap the outputs the notebook has stored |
| `closedform/...` | the docstring examples against their closed forms, computed as balls with Arb |

The criterion is always that the result is an **enclosure**: checks use
`overlaps`, never equality of balls. A value may move, as long as it still
encloses what it must.

`tests/test_dde.sage` covers `code/number_theoretic_dde_solutions.py`. It does not compare
with a recorded value: rather, it checks a known property.

| name | what it checks |
| --- | --- |
| `dde/closedform/...` | `ρ = 1` on `[0, 1]` and `1 − log u` on `[1, 2]`; `B = 0` below 1, `1` on `[1, 2]`, `1 + log(u−1)` on `[2, 3]`, at 53 and 150 bits |
| `dde/enclosure/...` | 30 bits meets 200 bits, at a spread of `u` up to 23/2 |
| `dde/continuity/...` | consecutive approximants agree at the step between them, up to the jump there: `+1` at 1 for `B`, `+1` at `α` and `−1` at 1 for `S_α` |
| `dde/jumps/...` | those two jumps of `S_α` measured over a window of `10^{-12}`, for `α = 2/5` and `1/2` |
| `dde/extension/...` | evaluating past the last step extends the grid and gives what a fresh instance gives, coefficient for coefficient |
| `dde/sanity/dickman-against-sage` | agreement with Sage's floating-point `dickman_rho` to `1e-12`; **not** an enclosure check, neither value bounding the other |
| `dde/errors/...`, `dde/arguments/...` | `α` outside `(0, 1)` or non-rational is refused; an `int`, a `Rational`, a real, a float and a `RealBall` argument all work |

`tests/test_timing.sage` holds micro-benchmarks of the parts the later steps
rewrite (the Taylor expansion and the LattE call on one balanced piece, the
elimination of the last coordinate, the Buchstab and Dickman solutions). They
are tagged `timing` and only run with `--only-timing`.

## Timing

```
sage tests/run_tests.sage --long --time after.json
sage tests/run_tests.sage --long --time after.json --compare before.json
```

`--time` writes the wall-clock time of every test that ran, with the commit,
the date and the machine; `--repeat N` (3 by default with `--time`) keeps the
median of N runs. `--compare` prints a before/after table with the ratios.
`--serial` replaces the library's pool of workers by an in-process map: a
serial time is only comparable with another serial time.

Timings mean nothing unless the machine is idle, and a before-and-after pair
has to be taken back to back, on the same machine, at the same setting.

## Reference values

`tests/reference_values.json` holds, for each example, the ball it returned at
its own precision, as exact midpoint and radius, with the commit that produced
it. Regenerate with

```
sage tests/make_references.sage               # all of them
sage tests/make_references.sage 'stadlmann/*' # only some, the rest kept
```

These are regression references, not certified values.

## Adding a test

Add a `tests/test_*.sage` file defining

```python
TESTS = {"group/name": (callable, ("quick",))}   # tags: quick, long, timing
KNOWN_ISSUES = {"group/name": "why it fails, and which step fixes it"}
```

A test raises `AssertionError` with a message that says what was expected and
what came out. The helpers of `tests/harness.sage` are already loaded:
`example_value(name, extra_bits=0)` caches values, `examples_namespace()`,
`h_namespace()` and `dde_namespace()` give access to the loaded scripts,
`in_directory(path)` runs a body elsewhere, and `ball_from_json` rebuilds a
stored ball.

Two things constrain how tests are written:

- the example scripts and `code/sieve_integral.py` are loaded into the globals
  of the test run, because the library's worker pool pickles its functions by
  name and only finds them in `__main__`. **Do not bind a name they use**, in
  particular the symbolic variables `t`, `u`, `v`, `w`, `z`, `x1` … `a7`, the
  Stadlmann script's dict `a`, or their `test_*` functions;
- `Computation of h/computation_of_h.sage` is loaded into a namespace of its
  own, since it would otherwise overwrite the library's
  `are_inequalities_compatible` with its own.
