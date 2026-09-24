# BCa Bootstrap Results

## Overview

This analysis used bias-corrected and accelerated (BCa) bootstrap intervals to quantify uncertainty in the two-layer functional mediation model:

```text
SNP exposure -> corpus callosum shape -> time to AD conversion
```

The analysis included 332 participants with baseline LMCI, 145 observed AD conversion events, and 187 censored observations. Four SNPs were evaluated in separate models. The corpus callosum square-root velocity function was treated as one multivariate functional mediator with two coordinate functions over 100 contour positions.

The adjustment set contained sex, standardized age, handedness, APOE4, genetic ancestry PC1 and PC2, and education years. The first-stage shape model also contained an intercept.

## Bootstrap Design

The BCa analysis used 5,000 subject-level pairs-bootstrap replicates with random seed `20260921`. The joint FPCA dimension was fixed at `K = 18` before resampling. This dimension explained 85.98% of shape variation in the original sample.

Each bootstrap replicate sampled 332 participants with replacement. The shape regression, 18-component joint FPCA basis, and Cox model were re-estimated after resampling. All four SNP models shared the same resampled participant indices within a replicate.

The acceleration correction used all 332 leave-one-subject-out jackknife fits. The complete two-layer model was re-estimated after removing each participant. Fixing `K = 18` prevented changes in the FPCA dimension from introducing a discrete model-selection step into the jackknife calculation.

## BCa Interval Construction

For each estimand, the bootstrap distribution was used to calculate the bias-correction term `z0`. The leave-one-subject-out estimates were used to calculate the acceleration term `a`. The nominal 2.5% and 97.5% probabilities were then transformed using `z0` and `a`, and the corresponding empirical bootstrap quantiles formed the BCa interval.

The acceleration estimates were small for the 12 principal direct, indirect, and total effects, ranging from -0.007 to 0.014. Most differences between the BCa and percentile intervals therefore came from the bias-correction term.

## Effect Definitions

The direct effect is the Cox coefficient for one additional coded allele:

```text
Direct HR = exp(tau).
```

The functional indirect effect combines both shape coordinates:

```text
Indirect HR = exp[integral alpha_1(s) beta_1(s) ds
                  + integral alpha_2(s) beta_2(s) ds].
```

The total effect is the product-method combination:

```text
Total HR = Direct HR x Indirect HR.
```

## Main Results

Hazard ratios are reported per one additional coded allele.

The prepared source object does not document the identity or orientation of the coded allele. The estimates should therefore be interpreted per additional coded allele until the allele metadata are verified.

| SNP | Direct HR (95% BCa CI) | Indirect HR (95% BCa CI) | Total HR (95% BCa CI) | Warning-free replicates |
|---|---:|---:|---:|---:|
| `rs929708` | 0.881 (0.628-1.296) | 1.070 (0.907-1.243) | 0.943 (0.660-1.383) | 4,998 |
| `rs11719939` | 1.050 (0.748-1.540) | 1.065 (0.885-1.246) | 1.118 (0.791-1.649) | 4,997 |
| `rs4639533` | 0.846 (0.597-1.212) | 1.131 (0.966-1.368) | 0.957 (0.665-1.361) | 4,998 |
| `rs2515029` | 0.981 (0.708-1.285) | 0.953 (0.834-1.102) | 0.934 (0.684-1.245) | 4,997 |

Every 95% BCa interval included HR = 1. The BCa analysis therefore did not identify a direct, shape-mediated, or total effect that was statistically distinguishable from the null at the 95% interval level.

`rs929708` and `rs4639533` had opposing direct and indirect point-estimate directions. `rs11719939` had direct and indirect HRs above 1, whereas `rs2515029` had both HRs below 1. These are descriptive pathway patterns and do not establish causal mediation.

`rs4639533` remained the closest to a positive shape-mediated association. Its indirect HR was 1.131 with a 95% BCa CI of 0.966-1.368.

## Indirect-Effect Corrections

| SNP | Bias correction `z0` | Acceleration `a` | Percentile 95% CI | BCa 95% CI |
|---|---:|---:|---:|---:|
| `rs929708` | -0.1800 | -0.0017 | 0.935-1.292 | 0.907-1.243 |
| `rs11719939` | -0.2279 | -0.0015 | 0.919-1.308 | 0.885-1.246 |
| `rs4639533` | -0.0005 | 0.0140 | 0.962-1.361 | 0.966-1.368 |
| `rs2515029` | 0.1473 | -0.0020 | 0.814-1.081 | 0.834-1.102 |

The small acceleration values indicate that the jackknife influence distribution contributed little additional skewness correction. The larger bias-correction terms for `rs929708`, `rs11719939`, and `rs2515029` shifted their BCa intervals relative to the uncorrected percentile intervals. These shifts did not change whether the intervals included 1.

## Computation and Diagnostics

All 5,000 bootstrap replicates completed without a replicate-level error. Ten individual Cox fits produced convergence warnings and were excluded only from the corresponding SNP summaries. The warning-free counts were 4,998 for `rs929708`, 4,997 for `rs11719939`, 4,998 for `rs4639533`, and 4,997 for `rs2515029`.

All 332 jackknife fits completed without an error or Cox warning. The fixed-`K` point estimates agreed with the original `K = 18` model estimates to numerical precision.

The reported intervals are model-specific 95% intervals. They were not adjusted for testing four SNPs. The analysis quantifies uncertainty in model-based associations and does not by itself establish causal direct or mediated effects.

## Output Files

| File | Content |
|---|---|
| [`bca_summary_B5000_seed20260921_K18.csv`](bca_summary_B5000_seed20260921_K18.csv) | Point estimates, bootstrap bias and standard errors, `z0`, acceleration, adjusted probabilities, BCa intervals, percentile intervals, and valid replicate counts |
| [`bca_bootstrap_estimates_B5000_seed20260921_K18.csv`](bca_bootstrap_estimates_B5000_seed20260921_K18.csv) | Replicate-level direct, coordinate-specific indirect, overall indirect, and total effects |
| [`bca_jackknife_estimates_B5000_seed20260921_K18.csv`](bca_jackknife_estimates_B5000_seed20260921_K18.csv) | Leave-one-subject-out estimates used to calculate acceleration |
| [`bca_point_estimates_B5000_seed20260921_K18.csv`](bca_point_estimates_B5000_seed20260921_K18.csv) | Original fixed-`K` point estimates |
| [`bca_bootstrap_warnings_B5000_seed20260921_K18.csv`](bca_bootstrap_warnings_B5000_seed20260921_K18.csv) | Cox convergence warnings from bootstrap replicates |
| [`bca_bootstrap_errors_B5000_seed20260921_K18.csv`](bca_bootstrap_errors_B5000_seed20260921_K18.csv) | Replicate-level errors; this file contains no failed replicates |
| [`bca_jackknife_warnings_B5000_seed20260921_K18.csv`](bca_jackknife_warnings_B5000_seed20260921_K18.csv) | Jackknife Cox warnings; this file contains no warnings |
| [`bca_jackknife_errors_B5000_seed20260921_K18.csv`](bca_jackknife_errors_B5000_seed20260921_K18.csv) | Jackknife errors; this file contains no failed fits |
| [`bca_results_B5000_seed20260921_K18.rds`](bca_results_B5000_seed20260921_K18.rds) | Complete R result object, including settings and diagnostics |

## Reproduction

Run the following command from the project root:

```bash
Rscript code/08_bca_bootstrap_four_snp_mediation.R 5000 20260921 18
```

The full mediation analysis overview is available in [`../README.md`](../README.md).
