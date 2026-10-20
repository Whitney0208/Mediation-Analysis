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

## Output Files

| File | Content |
|---|---|
| [`shape_volume_effect_summary.csv`](shape_volume_effect_summary.csv) | Point estimates of all effects, path coefficients, and path test statistics |
| [`cc_volume_mediator.csv`](cc_volume_mediator.csv) | Total CC volume, ICV, the ratio, and the standardized mediator for each participant |
| [`bootstrap_summary_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_summary_B10000_seed20260905_k20.csv) | Bootstrap standard errors, $z_0$, acceleration, BCa, percentile, and basic intervals |
| [`bootstrap_estimates_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_estimates_B10000_seed20260905_k20.csv) | Replicate-level estimates, including both path coefficients |
| [`jackknife_estimates_n332_k20.csv`](../shape_volume_mediation_bootstrap/jackknife_estimates_n332_k20.csv) | Leave-one-subject-out estimates for the BCa acceleration |
| [`bootstrap_warnings_B10000_seed20260905_k20.csv`](../shape_volume_mediation_bootstrap/bootstrap_warnings_B10000_seed20260905_k20.csv) | The two warning-flagged Cox fits |

## Reproduction

Run from the `version 2` directory. The last argument of script 23 is the number of CPU cores.

```bash
Rscript code/22_fit_shape_volume_mediation.R
Rscript code/23_bootstrap_shape_volume_mediation.R 10000 20260905 6
```

The full mediation analysis overview is available in [`../README.md`](../README.md).
