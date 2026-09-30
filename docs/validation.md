# Validation — 30 September 2026

## Software

Python 3.12; NumPy 2.5.3; SciPy 1.18.1; pytest 9.1.1; Matplotlib 3.11.2.
The package installed successfully using `pip install -e '.[dev,plots]'`.
Source and wheel distributions built successfully. The installed `twostage-be`
command produced valid JSON output.

## Release checks

All 46 tests passed on each supported Python version: **3.10, 3.11, 3.12,
3.13 and 3.14**, for a total of 230 successful test executions. Each version
also passed citation schema validation, wheel installation, command-line
execution, the analysis example, the five-method comparison example and plot
generation. The CI matrix includes exactly those five supported minor versions;
package metadata declares `>=3.10,<3.15`.

The MIT license, author metadata (Vangelis D. Karalis), citation, and keywords
were checked. Publication files contain no machine-specific paths or assistant
attribution. Generated caches, smoke results and installation metadata are
excluded from the source archive.

## Tests and executable examples

46 tests passed, including independent observation-level GLM comparisons for
one stage, two ordinary stages, and a second stage of two subjects.
Tests also cover interim branches, UL equality and rejection, planning
formulas, minimum even sample size in the Potvin engine, and seeded simulation.

Completed smoke grids:

| Preset | Scenarios | Trials per scenario | Purpose |
|---|---:|---:|---|
| comparison | 30 | 300 | All five methods run and export results |
| tables | 96 | 100 | All Table 1/2 scenarios and references are executable |
| figures | 187 | 5 | All caption-defined curves are executable |

Smoke grids check execution, not statistical reproduction. Plot generation
was tested and an output was visually inspected. The analysis example ran.
Smoke-run outputs are temporary diagnostics and are not included in the
published source package.

## Initial boundary comparison

True GMR=1.25, CVw=20%, N1=12, UL=62; 20,000 trials per method; seed=2026;
default legacy planning engine.

| Method | Estimated acceptance | 95% MC Wilson interval | Published percentage |
|---|---:|---:|---:|
| TSD-1 | 3.675% | 3.423%–3.945% | 3.80% |
| TSD-2 | 3.835% | 3.578%–4.110% | 3.95% |

The printed reference values lie within these intervals. This is preliminary
agreement for one scenario per method, not full validation of the article.
The study generator includes the inputs to reproduce both published tables
and compare percentages across all scenarios.

## Outstanding validation

- Full million-trial Table 1/2 runs and 100,000-trial figure runs.
- Same deterministic input data analysed in Python and a running MATLAB version.
- Large-simulation lower/upper boundary comparisons for all five methods.
- GitHub Actions execution after choosing and publishing the repository.
