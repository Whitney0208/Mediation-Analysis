# Wild Bootstrap Results

## Overview

This analysis used a subject-level exponential multiplier wild bootstrap to quantify uncertainty in the two-layer functional mediation model:

```text
SNP exposure -> corpus callosum shape -> time to AD conversion
```

The analysis retained the original 332 baseline LMCI participants, including 145 observed AD conversion events and 187 censored observations. Four SNPs were evaluated in separate models. The corpus callosum square-root velocity function was treated as one multivariate functional mediator with two coordinate functions over 100 contour positions.

The adjustment set contained sex, standardized age, handedness, APOE4, genetic ancestry PC1 and PC2, and education years. The first-stage shape model also contained an intercept.

## Bootstrap Design

The analysis generated 1,000 independent multiplier replicates with random seed `20260920`. For replicate `b`, participant `i` received an independent exponential weight:

```text
G_i^(b) ~ Exp(1),    w_i^(b) = G_i^(b) / mean(G^(b)).
```

The same subject weights were applied to the first-stage shape regression, joint weighted functional principal component analysis (FPCA), and weighted Cox partial likelihood. This joint weighting preserves the dependence between the SNP-to-shape and shape-to-hazard estimates. Subjects, genotype groups, observed survival outcomes, and Cox risk-set membership remained fixed.

The weighted FPCA basis was re-estimated in every replicate. The number of shape components was selected to explain at least 85% of weighted shape variation. The selected dimension was 16 in 188 replicates, 17 in 810 replicates, and 18 in 2 replicates.

This implementation is a perturbation or multiplier wild bootstrap. It is distinct from a residual sign-flipping bootstrap.

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

The reported total effect is the product-method combination:

```text
Total HR = Direct HR x Indirect HR.
```

The primary intervals are centered basic 95% multiplier intervals. Weighted-bootstrap percentile intervals are also stored in the summary output.

## Main Results

Hazard ratios are reported per one additional coded allele.

The prepared source object does not document the identity or orientation of the coded allele. The estimates should therefore be interpreted per additional coded allele until the allele metadata are verified.

| SNP | Direct HR (95% wild CI) | Indirect HR (95% wild CI) | Total HR (95% wild CI) | Warning-free replicates |
|---|---:|---:|---:|---:|
| `rs929708` | 0.881 (0.661-1.212) | 1.070 (0.909-1.213) | 0.943 (0.671-1.316) | 1,000 |
| `rs11719939` | 1.050 (0.788-1.476) | 1.065 (0.900-1.220) | 1.118 (0.794-1.570) | 1,000 |
| `rs4639533` | 0.846 (0.619-1.155) | 1.131 (0.956-1.300) | 0.957 (0.677-1.377) | 999 |
| `rs2515029` | 0.981 (0.770-1.261) | 0.953 (0.842-1.094) | 0.934 (0.714-1.231) | 998 |

All direct, indirect, and total 95% wild bootstrap intervals included HR = 1. The analysis therefore did not identify an effect that was statistically distinguishable from the null at the 95% interval level.

`rs929708` and `rs4639533` had direct HRs below 1 and indirect HRs above 1, indicating opposing pathway directions in the point estimates. `rs11719939` had direct and indirect HRs above 1. `rs2515029` had direct and indirect HRs below 1. These patterns describe point-estimate directions and do not establish causal mediation.

Among the four SNPs, `rs4639533` had the largest positive indirect point estimate. Its indirect HR was 1.131 with a 95% wild CI of 0.956-1.300, which still included 1.

## Comparison With the Pairs Bootstrap

Wild bootstrap standard errors were smaller than the corresponding subject-level pairs-bootstrap standard errors for all 12 principal effects. For the indirect effects, the standard-error reductions were approximately 3% for `rs929708`, 4% for `rs11719939`, 9% for `rs4639533`, and 6% for `rs2515029`.

The modest reduction in uncertainty did not change the inferential conclusion. Both methods produced 95% intervals containing HR = 1 for every direct, indirect, and total effect.

## Computation and Diagnostics

All 1,000 multiplier replicates completed, and no replicate failed. Three individual Cox fits produced convergence warnings. Two warnings occurred for `rs2515029`, and one occurred for `rs4639533`. These fits were omitted only from the corresponding SNP summaries. Results for `rs929708` and `rs11719939` used all 1,000 replicates.

The exponential multiplier weights had a median effective sample size of approximately 168. This is expected because weighting changes the information contribution of participants while retaining all observed records.

The reported intervals are model-specific 95% intervals. They were not adjusted for testing four SNPs. The analysis supports uncertainty assessment for model-based associations and does not by itself establish causal direct or mediated effects.

## Output Files

| File | Content |
|---|---|
| [`wild_bootstrap_summary_B1000_seed20260920.csv`](wild_bootstrap_summary_B1000_seed20260920.csv) | Point estimates, wild standard errors, primary basic intervals, secondary percentile intervals, and valid replicate counts |
| [`wild_bootstrap_estimates_B1000_seed20260920.csv`](wild_bootstrap_estimates_B1000_seed20260920.csv) | Replicate-level direct, coordinate-specific indirect, overall indirect, and total effects |
| [`wild_bootstrap_weight_diagnostics_B1000_seed20260920.csv`](wild_bootstrap_weight_diagnostics_B1000_seed20260920.csv) | Effective sample size and multiplier range for every replicate |
| [`wild_bootstrap_warnings_B1000_seed20260920.csv`](wild_bootstrap_warnings_B1000_seed20260920.csv) | Cox convergence warnings |
| [`wild_bootstrap_errors_B1000_seed20260920.csv`](wild_bootstrap_errors_B1000_seed20260920.csv) | Replicate-level errors; this file contains no failed replicates |
| [`wild_bootstrap_results_B1000_seed20260920.rds`](wild_bootstrap_results_B1000_seed20260920.rds) | Complete R result object, including settings and diagnostics |

## Reproduction

Run the following command from the project root:

```bash
Rscript code/07_wild_bootstrap_four_snp_mediation.R 1000 20260920
```

The full mediation analysis overview is available in [`../README.md`](../README.md).
