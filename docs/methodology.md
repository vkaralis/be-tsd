# Methodology and source mapping

## Scope

One endpoint, two treatments, two periods, two sequences, complete observations,
and equal sequence allocation in each stage. Both stages share a treatment
effect and within-subject variance. Dropouts, unequal sequence sizes, replicate
designs, scaled BE and multi-endpoint correlations are outside this implementation.

Input arrays have shape `(n/2, 2)` with TR and RT sequence columns. Each entry is
the subject's **ln(Test) - ln(Reference)**. Do not pass raw AUC/Cmax values.
The generator simulates differences directly with mean ln(true GMR) and
variance `2 * log(1 + CVw**2)`. This has the same difference distribution as
the independent log-PK observations in the supplied MATLAB scripts.

## Balanced ANOVA reduction

Subject intercepts cancel in paired differences. Sequence, stage and nested
subject intercepts are absorbed by the fixed subject effects in the full model.
The difference regression retains a common treatment mean and one sequence
contrast per stage, representing Period(Stage). In a balanced design, its
treatment estimate is the overall mean difference. Divide its residual sum of
squares by two to obtain the observation-scale ANOVA residual sum of squares.

Residual df is `N1 - 2` after stage 1 and `N1 + N2 - 3` after combining stages.
The combined residual includes between-stage treatment-mean heterogeneity;
simply pooling the two within-stage residual variances would be incorrect.
The standard error is `sqrt(2 * MSE / N)`. The log-scale CI is
`log_GMR +/- t(1-alpha, df) * SE`; BE requires both ends within ln(0.80)-ln(1.25).

## Planning engines

`potvin` uses the central-t CDF difference in Potvin section 2.3. Required N is
the smallest even total size with estimated power >= target, using planning
df=N-2. The actual combined test uses N-3. This is the paper's planning
approximation; it is not an exact unconditional TOST power calculation.

`legacy` reproduces the nearest-boundary formula found in the supplied MATLAB:

```text
d = min(log(GMR / 0.80), log(1.25 / GMR))
q = t_quantile(1-alpha, N-2)
power = t_cdf(-q + sqrt((N-q/2) * d**2 / (2*MSE)), N-2)
```

Sample size starts from the normal-quantile approximation
`2*MSE*(za+zb)**2/d**2 + za/2`, then iterates using t quantiles and df=round(N)-2.
The supplied code's one-sided tolerance (`new_N - old_N < 0.5`) is retained.
At GMR exactly 1, the beta quantile is based on beta/2. The final size is rounded
as MATLAB round-half-up, then increased to even. A minimum total N=4 is enforced.
The loop exits at convergence and reports failure if 100 iterations do not
converge. Legacy estimated power and sample-size formulas are distinct
approximations, so the resulting N need not be the minimum attaining 80%
under the legacy power formula. Defaults preserve the authors' MATLAB planning
formula for TSD while using the published Potvin planning formula for B/C/D.

`power_engine` can be overridden for sensitivity analyses. Changing it changes
the design's operating characteristics; simulation output always records it.

## Source mapping and corrections

| Source | Role |
|---|---|
| Karalis (2013), sections 2.1.1-2.1.3 | TSD-1/2 branches, GMR gate, alpha and UL |
| `7_TS_BCD_gmr_EMA/B_gmr_EMA/mainTSD.m` | TSD-2 BE-first branch and fixed-effect analysis |
| `1_V_method/0_GUI/coreV.m`, `6_TS_CD_gmr/C_gmr/coreV.m` | TSD-1 power-first branch and planning expressions |
| `5_TS_B_gmr/N2.m` | Julious-style legacy sample-size iteration |
| `2_Potvin_method/0_GUI/core.m` | Fixed planning GMR in C/D; historical implementations differ |
| Potvin (2008), sections 2.2-2.3 | B/C/D rules and planning approximation |
| Montague (2012) | Default planning GMR=.90 for D |

These paths are relative to the original parent workspace and are not required
to run the Python package. No MATLAB executable was run during initial validation.

Explicit choices where historical scripts differ:

- **UL equality is allowed**, following the article's N <= UL. Some scripts use
  `nadd >= naddmax`, incorrectly excluding equality relative to that definition.
- UL is a rejection boundary, not a recruitment cap: if required N exceeds UL,
  stage 2 is not generated and the study fails with total recruitment N1.
- A continuing trial recruits at least two new subjects. N1 >= UL permits
  assessment at stage 1 but prevents any second-stage recruitment.
- The article's GMR gate includes endpoints. Exact endpoints have no finite
  sample-size solution for positive residual variance under these planning rules.
- Potvin B/C/D have neither the observed-GMR gate nor a mandatory finite UL.
  The supplied historical C/D scripts include cap-related modifications.
- Every actual analysis uses fixed subject effects and Period(Stage); old scripts
  using random-subject settings are not blindly copied.
- The figure preset uses N1=24/48 from the captions, rather than the conflicting
  N1=62/72/150 sentence in section 2.3.

## Published reference values

`article_tables.csv` transcribes Tables 1 and 2 as printed (percent units).
The tables preset assumes true GMR=1.25 for boundary type I error, consistent
with the MATLAB comment. Repeated entries in the source table are preserved.
Reference values are comparisons, not hard-coded expectations that every
modern implementation must reproduce exactly.

Seeds use NumPy Generator; study grids derive independent scenario seeds with
SeedSequence.spawn. Use the same software versions as well as seed and settings
when requiring identical numerical output. CSV and JSON results record settings.

The initial validation covers algebra and decisions. Full-grid statistical
validation, lower-boundary symmetry checks at large simulation counts and
deterministic same-data comparison with MATLAB remain future validation work.
