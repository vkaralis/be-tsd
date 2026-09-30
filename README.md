# Two-stage bioequivalence designs in Python

**Author:** Vangelis D. Karalis  
**License:** [MIT](LICENSE)  
**Tags:** Bioequivalence, Two-stage designs

Python implementation of Karalis **TSD-1 / TSD-2** and **Potvin B / C / D**
for complete, balanced 2 x 2 crossover studies with one PK endpoint.
Includes interim decisions, sample-size re-estimation, combined-stage analysis,
seeded Monte Carlo simulations and article scenario grids.

## Article citation

This implementation of TSD-1/TSD-2 is based on the following article. Please cite it when using these methods:

> Karalis, V. (2013). The role of the upper sample size limit in two-stage bioequivalence designs. *International Journal of Pharmaceutics*, **456**(1), 87–94. [https://doi.org/10.1016/j.ijpharm.2013.08.013](https://doi.org/10.1016/j.ijpharm.2013.08.013)

## Install

Supported Python versions: **3.10, 3.11, 3.12, 3.13 and 3.14**.
The CI matrix covers every supported version. From this directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev,plots]'
python -m pytest
```

## Quick start

```bash
twostage-be --method TSD-1 --n1 24 --gmr 0.95 --cv 0.20 \
  --upper-limit 150 --trials 10000 --seed 2026 --output results/tsd1.json
twostage-be --method B --n1 24 --gmr 1.25 --cv 0.20 --trials 10000
twostage-be --method C --n1 24 --gmr 0.95 --cv 0.40 --trials 10000
twostage-be --method D --n1 24 --gmr 0.90 --cv 0.40 --trials 10000
```

`--gmr` is the true population ratio for generating data.
`--assumed-gmr` is the fixed planning ratio for Potvin methods only.
CV and probabilities are fractions: `0.20` means 20%.

```python
from twostage_be import Design, simulate

result = simulate(Design("TSD-2", upper_limit=150),
                  n1=24, gmr=0.95, cv=0.20, trials=1000, seed=2026)
print(result["acceptance_probability"], result["mean_total_n"])
```

## Methods

| Method | Planning GMR | Initial decision | Adjusted alpha | Default planning engine |
|---|---|---|---|---|
| TSD-1 | Observed | Power at .05; if >=80%, test at .05 and stop; otherwise test at .028 | .028 | legacy |
| TSD-2 | Observed | Test at .0294; if unsuccessful, check power at .0294 | .0294 | legacy |
| B | Fixed .95 | Test at .0294; if unsuccessful, check power at .0294 | .0294 | potvin |
| C | Fixed .95 | Power at .05; if >=80%, test at .05 and stop; otherwise test at .0294 | .0294 | potvin |
| D | Fixed .90 | Same decision structure as C, adjusted alpha .028 | .028 | potvin |

TSD-1/2 first reject an observed GMR outside 0.80-1.25 and require a predefined
total upper limit. Potvin methods have no observed-GMR gate and no upper limit
by default. Setting a finite UL for Potvin produces a capped variant.
All methods recruit at least two additional subjects if they continue.

## Article grids and plots

```bash
python examples/run_study.py --preset comparison --trials 1000
python examples/plot_study.py results/study.csv
python examples/analyse_trial.py

# Full grids are computationally expensive; start with a small --trials value.
python examples/run_study.py --preset tables --trials 1000000 --output results/tables.csv
python examples/run_study.py --preset figures --trials 100000 --output results/figures.csv
python examples/plot_study.py results/figures.csv
```

The tables preset covers the 96 scenarios in Tables 1/2. It includes published
percentages and differences in percentage points. The figures preset follows
the captions for Figures 2/3/5 and adds the N1=72, UL=72 curve from Figure 4.
The paragraph in section 2.3 disagrees with the figure captions about N1;
these choices are explicit in [methodology](docs/methodology.md).

The acceptance probability at true GMR 0.80 or 1.25 estimates boundary type I
error. Each scenario reports a Monte Carlo standard error, a 95% Wilson interval,
stage-specific acceptance, second-stage frequency, mean total N and stop reasons.
The interval describes simulation uncertainty, not the BE confidence interval.

## Validation status

Tests check the paired-difference analysis against an independently constructed
observation-level fixed-subject GLM, including a second stage of only two
subjects. They also check decision branches, sample-size minimality for the
Potvin engine, UL equality, numerical formulas, invalid inputs and reproducibility.

This is an initial research implementation. Full million-trial reproduction and
same-input numerical comparison with a running MATLAB implementation have not
yet been completed. Passing tests does not establish type I error control for
every parameter combination. See [methodology](docs/methodology.md) for the
legacy formula and documented differences from the supplied MATLAB versions.
The completed checks and initial numerical comparison are recorded in
[validation](docs/validation.md).

## References

Please cite the following article when using the TSD-1/TSD-2 implementation:

- Karalis, V. (2013). The role of the upper sample size limit in two-stage bioequivalence designs. *International Journal of Pharmaceutics*, **456**(1), 87–94. [https://doi.org/10.1016/j.ijpharm.2013.08.013](https://doi.org/10.1016/j.ijpharm.2013.08.013)

GitHub's **Cite this repository** feature uses [CITATION.cff](CITATION.cff),
which lists this article as the preferred citation.

For **Potvin B/C**, cite Potvin et al.; for **D with assumed GMR=0.90**,
also cite Montague et al. The package implements all three methods alongside
TSD-1/TSD-2:

- Potvin D. et al. (2008). *Sequential design approaches for bioequivalence
  studies with crossover designs.* Pharmaceutical Statistics 7, 245-262.
  [doi:10.1002/pst.294](https://doi.org/10.1002/pst.294)
- Montague T.H. et al. (2012). *Additional results for Sequential design
  approaches for bioequivalence studies with crossover designs.*
  Pharmaceutical Statistics 11, 8-13.
  [doi:10.1002/pst.483](https://doi.org/10.1002/pst.483)

## License

Copyright (c) 2026 Vangelis D. Karalis. Released under the [MIT License](LICENSE).
