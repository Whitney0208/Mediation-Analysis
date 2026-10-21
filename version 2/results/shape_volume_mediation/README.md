# Shape and Corpus Callosum Volume as Joint Mediators

## Overview

This supplementary analysis adds total corpus callosum (CC) volume to the main penalized functional Cox model, so that shape and size enter the same mediation model:

```text
SNP exposure -> CC shape (SRVF)  -> time to AD conversion
SNP exposure -> CC volume / ICV  -> time to AD conversion
```

It addresses whether the shape-mediated pathway simply reflects CC size. The SRVF shape representation removes scale, so CC volume carries the size information that the shape mediator does not.

The sample, adjustment set, and shape model are the same as in the [main analysis](../penalized_functional_cox/README.md): 332 participants with baseline LMCI, 145 AD conversion events, and adjustment for sex, standardized age, handedness, APOE4, genetic ancestry PC1 and PC2, and education years.

## Volume mediator

The five FreeSurfer CC subregion volumes in `data/SPLS-Coxdata/CCinfo.csv` (anterior, mid-anterior, central, mid-posterior, posterior) were summed and divided by intracranial volume (ICV) from the same file:

$$V_i = \frac{\sum_{r=1}^{5}\text{volume}_{ir}}{\text{ICV}_i},$$

then standardized to mean 0 and SD 1 in the 332 participants. All participants matched by RID and had complete values. The standardization is a fixed linear transformation and does not affect the indirect effects.

## Model

### First layer

$$M_j(s) = \alpha_j(s)\,X + \lambda_j(s)^\top W + \varepsilon_j(s), \qquad V = a\,X + c^\top W + e$$

The shape layer is unchanged: pointwise least squares followed by GCV smoothing splines. The volume layer is an ordinary least-squares regression with an intercept.

### Second layer

$$h(t) = h_0(t)\exp\Big[\tau X + \kappa^\top W + \sum_{j=1}^{2}\int_0^1 M_j(s)\,\beta_j(s)\,ds + b\,V\Big]$$

Each $\beta_j(s)$ uses 20 cubic P-spline basis functions with a second-order difference penalty and REML smoothing, exactly as in the main model. The volume coefficient $b$ is unpenalized.

### Effects

| Effect | Definition |
|---|---|
| Direct | $\tau$ |
| Indirect through shape | $\sum_j \int \alpha_j(s)\beta_j(s)\,ds$ |
| Indirect through volume | $a \times b$ |
| Total indirect | Shape + volume |
| Total | Direct + total indirect |

Separating the shape and volume pathways assumes that neither mediator causally affects the other. The total indirect effect does not require this assumption.

## Inference

Pairs bootstrap with 10,000 replicates and seed `20260905`, the same subject indices as the main pairs bootstrap. Both layers were refitted in every replicate, with REML reselection of the smoothing parameters. BCa acceleration used all 332 leave-one-subject-out fits. No replicate failed. Two replicate-level Cox fits produced warnings and were excluded from the corresponding SNP summaries, so some effects are based on 9,999 replicates. No jackknife fit produced a warning.

## Results

Hazard ratios are per one additional coded allele.

### Indirect effects (95% BCa CI)

| SNP | Through shape | Through volume | Total indirect | Main model, shape only |
|---|---:|---:|---:|---:|
| `rs929708` | 1.092 (0.981-1.223) | 0.999 (0.967-1.013) | 1.091 (0.976-1.221) | 1.095 (0.985-1.224) |
| `rs11719939` | 1.088 (0.972-1.213) | 0.999 (0.968-1.016) | 1.087 (0.967-1.215) | 1.091 (0.975-1.219) |
| `rs4639533` | 1.131 (0.993-1.317) | 1.003 (0.989-1.051) | 1.133 (0.994-1.322) | 1.135 (0.997-1.322) |
| `rs2515029` | 0.952 (0.856-1.048) | 0.996 (0.951-1.015) | 0.948 (0.853-1.048) | 0.950 (0.857-1.045) |

Under the percentile interval, the `rs4639533` shape and total indirect effects were 1.005-1.340 and 1.002-1.340. Every basic interval included 1.

### Direct and total effects (95% BCa CI)

| SNP | Direct HR | Total HR |
|---|---:|---:|
| `rs929708` | 0.840 (0.620-1.165) | 0.916 (0.665-1.280) |
| `rs11719939` | 0.999 (0.741-1.388) | 1.085 (0.796-1.518) |
| `rs4639533` | 0.864 (0.643-1.180) | 0.979 (0.706-1.365) |
| `rs2515029` | 1.008 (0.772-1.312) | 0.956 (0.730-1.260) |

### Volume pathway components

| SNP | SNP to volume: t | Volume to log hazard: z | Volume to log hazard: p |
|---|---:|---:|---:|
| `rs929708` | 0.36 | -0.37 | 0.71 |
| `rs11719939` | 0.25 | -0.36 | 0.72 |
| `rs4639533` | -0.79 | -0.37 | 0.71 |
| `rs2515029` | 1.32 | -0.36 | 0.72 |

## Interpretation

1. CC volume did not mediate the SNP effects. Neither path of the volume pathway was associated, and the volume indirect HRs were 0.996 to 1.003.
2. Adding volume left the shape results essentially unchanged. Shape indirect HRs moved by about 0.003, and standard errors increased by 1% to 3%.
3. The mediation signal therefore reflects CC shape rather than CC size. This agrees with [S-GMAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC12311987/), which found no mediation through the five CC subregion volumes but did find shape-mediated pathways.
4. As in the main analysis, `rs4639533` shows a borderline shape-mediated effect that is not robust to the choice of interval.

For `rs4639533`, the BCa interval for the volume pathway is shifted relative to the percentile interval ($z_0 = 0.43$). The product of two near-zero coefficients has a sharply peaked bootstrap distribution, for which BCa corrections are less reliable. All volume-pathway intervals include 1 under every interval type.

## Wild Bootstrap

For direct comparison with the wild bootstrap of the main model, the same model was also analyzed with the multiplier wild bootstrap: 10,000 replicates with seed `20260920` (the same Exp(1) subject weights as the main-model wild bootstrap), weights applied to the shape layer, the volume regression, and the penalized Cox layer, and smoothing parameters held at their full-data REML values. BCa acceleration used 332 leave-one-subject-out fits with the same fixed smoothing parameters. No replicate failed and no fit produced a warning.

| SNP | Effect | Basic | Percentile | BCa | Pairs BCa (for reference) |
|---|---|---:|---:|---:|---:|
| `rs929708` | Through shape | 0.980-1.192 | **1.000-1.217** | **1.002-1.219** | 0.981-1.223 |
| | Total indirect | 0.977-1.196 | 0.995-1.217 | 0.996-1.219 | 0.976-1.221 |
| `rs11719939` | Through shape | 0.980-1.187 | 0.997-1.207 | 1.000-1.212 | 0.972-1.213 |
| | Total indirect | 0.977-1.192 | 0.991-1.209 | 0.994-1.213 | 0.967-1.215 |
| `rs4639533` | Through shape | 0.986-1.262 | **1.013-1.296** | **1.016-1.301** | 0.993-1.317 |
| | Total indirect | 0.988-1.272 | **1.010-1.300** | **1.015-1.307** | 0.994-1.322 |
| `rs2515029` | Through shape | 0.878-1.052 | 0.861-1.032 | 0.863-1.034 | 0.856-1.048 |
| | Total indirect | 0.872-1.049 | 0.857-1.031 | 0.859-1.032 | 0.853-1.048 |

Values are HRs with 95% intervals; bold intervals exclude 1. Exact lower bounds that round to 1.000: `rs929708` shape percentile 1.0003 (excludes 1), `rs11719939` shape BCa 0.9997 (includes 1). The volume pathway remained null under every interval (indirect HRs 0.996 to 1.003), and every direct- and total-effect interval included 1.

The pattern matches the shape-only model: under the wild percentile and BCa intervals, the shape-mediated effects of `rs4639533` and `rs929708` exclude 1, whereas the basic intervals and the pairs-bootstrap BCa intervals include 1. Because the wild bootstrap does not propagate smoothing-parameter uncertainty, its standard errors were 12% to 15% smaller than the pairs-bootstrap standard errors for the shape and total indirect effects, so the pairs bootstrap remains the primary analysis.

## Sensitivity: ICV as a Covariate

The main version divides total CC volume by ICV. In this sensitivity analysis, the mediator is the raw total CC volume (standardized), and standardized ICV is added to the adjustment set of every equation in both layers. The bootstrap design is identical (10,000 pairs replicates, seed `20260905`, BCa from 332 jackknife fits; two warning-flagged fits excluded).

| SNP | Through shape: ratio | Through shape: ICV covariate | Total indirect: ratio | Total indirect: ICV covariate |
|---|---:|---:|---:|---:|
| `rs929708` | 1.092 (0.981-1.223) | 1.084 (0.972-1.215) | 1.091 (0.976-1.221) | 1.079 (0.962-1.213) |
| `rs11719939` | 1.088 (0.972-1.213) | 1.079 (0.962-1.202) | 1.087 (0.967-1.215) | 1.073 (0.953-1.201) |
| `rs4639533` | 1.131 (0.993-1.317) | 1.134 (1.000-1.311) | 1.133 (0.994-1.322) | 1.136 (0.997-1.315) |
| `rs2515029` | 0.952 (0.856-1.048) | 0.955 (0.857-1.052) | 0.948 (0.853-1.048) | 0.951 (0.853-1.053) |

Values are HRs with 95% BCa intervals. The volume pathway remained null (indirect HRs 0.995 to 1.002), and direct and total effects were essentially unchanged. The shape indirect effect of `rs4639533` had a BCa lower bound of 1.0004 in this version, compared with 0.993 in the ratio version, so the borderline status of this SNP does not depend on how ICV is handled.

## Volume-Only Mediation (Comparison with S-GMAS)

To compare directly with the volumetric analysis of S-GMAS, a separate model used CC volumes as the only mediators, without the shape term. Each of the five subregion volumes was divided by ICV and standardized. The second layer is an ordinary Cox model, and the bootstrap design matches the main analysis (10,000 pairs replicates, seed `20260905`, BCa from 332 jackknife fits; 8 of 40,000 replicate-level fits produced Cox warnings and were excluded).

- **Joint model.** All five subregions enter together: $V_r = a_r X + c_r^\top W + e_r$ and $h(t) = h_0(t)\exp[\tau X + \kappa^\top W + \sum_r b_r V_r]$, with indirect effect $\sum_r a_r b_r$.
- **Single-region models.** One model per subregion, with indirect effect $a_r b_r$.

| SNP | Joint indirect | Anterior | Mid-anterior | Central | Mid-posterior | Posterior |
|---|---:|---:|---:|---:|---:|---:|
| `rs929708` | 1.016 (0.927-1.105) | 0.981 (0.925-1.004) | 1.003 (0.963-1.051) | 1.005 (0.938-1.064) | 1.017 (0.971-1.081) | 0.998 (0.968-1.010) |
| `rs11719939` | 0.982 (0.884-1.065) | 0.984 (0.927-1.010) | 0.999 (0.950-1.045) | 0.998 (0.918-1.060) | 0.998 (0.940-1.049) | 1.002 (0.988-1.036) |
| `rs4639533` | 1.056 (0.972-1.180) | 0.999 (0.953-1.039) | 1.020 (0.982-1.090) | 1.052 (0.995-1.148) | 1.027 (0.985-1.104) | 1.004 (0.992-1.043) |
| `rs2515029` | 0.966 (0.898-1.037) | 0.995 (0.957-1.020) | **0.965 (0.911-0.997)** | 0.963 (0.903-1.004) | 0.970 (0.914-1.003) | 0.994 (0.954-1.008) |

Values are indirect HRs with 95% BCa intervals. No joint volumetric indirect effect excluded 1. One of the 20 single-region tests (`rs2515029`, mid-anterior; bold) excluded 1, which is about the number expected by chance at the 5% level, and it would not survive a Bonferroni correction for 20 tests (threshold 0.0025). This agrees with S-GMAS, which found no mediation through either the total or the individual CC subregion volumes.

## Output Files

| File | Content |
|---|---|
| [`shape_volume_effect_summary.csv`](shape_volume_effect_summary.csv) | Point estimates of all effects, path coefficients, and path test statistics |
| [`cc_volume_mediator.csv`](cc_volume_mediator.csv) | Total CC volume, ICV, the ratio, and the standardized mediator for each participant |
| [`bootstrap_summary_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_summary_B10000_seed20260905_k20.csv) | Bootstrap standard errors, $z_0$, acceleration, BCa, percentile, and basic intervals |
| [`bootstrap_estimates_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_estimates_B10000_seed20260905_k20.csv) | Replicate-level estimates, including both path coefficients |
| [`jackknife_estimates_n332_k20.csv`](../shape_volume_mediation_bootstrap/jackknife_estimates_n332_k20.csv) | Leave-one-subject-out estimates for the BCa acceleration |
| [`bootstrap_warnings_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_warnings_B10000_seed20260905_k20.csv) | The two warning-flagged Cox fits |
| [`wild_bootstrap_summary_B10000_seed20260920_k20.csv`](../shape_volume_wild_bootstrap/wild_bootstrap_summary_B10000_seed20260920_k20.csv) | Wild bootstrap: standard errors, $z_0$, acceleration, basic, percentile, and BCa intervals |
| [`wild_bootstrap_estimates_B10000_seed20260920_k20.csv`](../shape_volume_wild_bootstrap/wild_bootstrap_estimates_B10000_seed20260920_k20.csv) | Replicate-level wild bootstrap estimates |
| [`fixed_smoothing_jackknife_n332_k20.csv`](../shape_volume_wild_bootstrap/fixed_smoothing_jackknife_n332_k20.csv) | Leave-one-subject-out estimates for the wild BCa acceleration |
| [`bootstrap_summary_B10000_seed20260905_k20.csv`](../shape_volume_icv_covariate_bootstrap/bootstrap_summary_B10000_seed20260905_k20.csv) | ICV-covariate sensitivity: point estimates and intervals |
| [`bootstrap_summary_B10000_seed20260905.csv`](../volume_only_mediation/bootstrap_summary_B10000_seed20260905.csv) | Volume-only mediation: joint and single-region indirect effects with BCa and percentile intervals |
| [`volume_only_point_estimates.csv`](../volume_only_mediation/volume_only_point_estimates.csv) | Volume-only point estimates |

## Reproduction

Run from the `version 2` directory. The last argument of scripts 23, 26, and 27 is the number of CPU cores.

```bash
Rscript code/22_fit_shape_volume_mediation.R
Rscript code/23_bootstrap_shape_volume_mediation.R 10000 20260905 6
Rscript code/28_wild_bootstrap_shape_volume.R 10000 20260920
Rscript code/27_shape_volume_icv_covariate.R 10000 20260905 6
Rscript code/26_volume_only_mediation.R 10000 20260905 6
```

The full mediation analysis overview is available in [`../README.md`](../README.md).
