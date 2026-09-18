# Examples

## Content

This folder containes SageMath Jupyter notebooks providing the details for the main examples described in the preprint *Computing sieve integrals using LattE, and the density of integers with a localized divisor* ([arXiv](https://arxiv.org/abs/2606.30428), [HAL](https://hal.science/hal-05673627)).

## Usage

* [Launch Sage](https://doc.sagemath.org/html/en/installation/launching.html) with the jupyter notebook:
```
sage -n jupyter
```
* Open the example file.
* Execute cells from top to bottom.
* Feel free to vary the ``precision``, ``facteur`` and ``verbose`` parameters.

## Scripted versions of the same tests

The same integrals are also available as plain SageMath scripts, one per
notebook, so that they can be run in batch, timed and profiled without going
through Jupyter. They are ``.sage`` files, preparsed like the notebook cells,
so that ``1/2`` is a rational and a floating point literal is a ``RealNumber``
exactly as there. As with the notebooks, run them **from this directory**:

```
cd Examples

# compute every example of one notebook
sage tests_chen.sage

# list the examples, then time a selection (shell globs on the names)
sage run_examples.sage --list
sage run_examples.sage 'chen/*' 'merikoski/F[123]'

# profile one example
sage run_examples.sage --profile --limit 30 maynard/I1
```

``run_examples.sage`` collects every example into ``ALL_EXAMPLES``, a
dictionary mapping a name such as ``merikoski/F4`` to a callable taking
``precision`` and ``verbose`` keyword arguments; each example's default
precision is the one its source notebook uses. From a Sage session started
here:

```python
sage: load("run_examples.sage")
sage: run("stadlmann/1*")
sage: profile("maynard/I1", precision=25)
```

Note that ``sieve_integral`` hands the pieces of a balanced polytope to a
``multiprocessing.Pool`` once there are more than twenty of them, and cProfile
only sees the parent process. ``--profile`` therefore runs the computation in a
single process by default (``--parallel`` keeps the pool, ``--serial`` forces
the single process for a plain timing run); wall clock timings taken the two
ways are not comparable.

Two examples are left out of a bare run over everything and have to be named
explicitly: ``stadlmann/40`` and ``stadlmann/41``, the six-dimensional
integrals on page 67 of Stadlmann's paper, which take substantially longer than
the rest.

## Detailed description of content

### ``Tests-Chen.ipynb``

Tests on two integrals from J. R. Chen's 1978 work on primes with a twin almost-prime (*Sci. Sin.* **21**, 421-430 (1978)).

### ``Tests-Ford-Maynard.ipynb``

Tests on five integrals from K. Ford and J. Maynard, [On the theory of prime producing sieves](https://arxiv.org/abs/2407.14368).

### ``Tests-Maynard.ipynb``

Tests on the integrals $`I_1, \dotsc I_9`$ from J. Maynard, [Primes with restricted digits](https://doi.org/10.1007/s00222-019-00865-6).

### ``Tests-Merikoski.ipynb``

Tests on the integrals $`F_1, \dotsc, F_6`$ from J. Merikoski, [On the largest prime factor of n^2 + 1](https://doi.org/10.4171/JEMS/1216).

### ``Tests-Stadlmann.ipynb``

Tests on several integrals from J. Stadlmann, [On the mean square gap between primes](https://arxiv.org/abs/2212.10867).

## Help

Please contact [@sarydrappeau](https://github.com/sarydrappeau) for help using these notebooks.

