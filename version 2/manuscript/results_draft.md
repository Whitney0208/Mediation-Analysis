<!--
Draft of the Results section. Items marked [TODO] must be resolved before submission:
  1. Primary interval: this draft uses the wild-bootstrap BCa interval; percentile intervals are given in Table 3.
  2. Time unit of the survival scale (assumed to be days).
  3. Anatomical labels of the contour ends (s = 0 and s = 0.5).
  4. Description of alpha(s) from Figure 1.
All effects are hazard ratios (HRs) per additional coded allele.
-->

# Results

## Study sample

The analysis included 332 participants with LMCI at baseline (Table 1). During follow-up, 145 (43.7%) converted to AD and 187 (56.3%) were censored. Median follow-up, estimated by the reverse Kaplan–Meier method, was 1099 days [TODO: confirm unit], and the Kaplan–Meier median time to conversion was 1121 days (95% CI 1094 to not reached). The cumulative incidence of conversion was 10.6% at 1 year, 34.7% at 2 years, and 44.8% at 3 years.

Mean age was 74.9 years (SD 7.1), 64.2% of participants were male, and all were White. Participants who converted were more often APOE ε4 carriers than those who did not (71.7% vs 43.9%) and had poorer baseline cognition (mean MMSE 26.7 vs 27.4; mean ADAS-Cog 13 21.2 vs 16.5). Age, sex, education, and corpus callosum (CC) volume were similar between the two groups.

Genotype frequencies of all four SNPs were consistent with Hardy–Weinberg equilibrium (exact p = 0.43–0.77; Table S1). rs929708 and rs11719939 were in strong linkage disequilibrium (r² = 0.75), whereas all other pairs were nearly uncorrelated (r² ≤ 0.008). Results for rs929708 and rs11719939 therefore reflect largely the same genetic signal. The coded allele of rs4639533 was the major allele (frequency 0.795).

**Table 1. Baseline characteristics.**

| Characteristic | Overall (n = 332) | Converted (n = 145) | Did not convert (n = 187) |
|---|---|---|---|
| Age, years, mean (SD) | 74.9 (7.1) | 74.7 (6.7) | 75.0 (7.5) |
| Male, n (%) | 213 (64.2) | 88 (60.7) | 125 (66.8) |
| Education, years, mean (SD) | 15.7 (3.0) | 15.6 (3.0) | 15.7 (3.0) |
| Right-handed, n (%) | 303 (91.3) | 133 (91.7) | 170 (90.9) |
| White, n (%) | 332 (100.0) | 145 (100.0) | 187 (100.0) |
| Not Hispanic or Latino, n (%) | 324 (97.6) | 143 (98.6) | 181 (96.8) |
| APOE ε4 alleles, n (%) | | | |
| &nbsp;&nbsp;0 | 146 (44.0) | 41 (28.3) | 105 (56.1) |
| &nbsp;&nbsp;1 | 143 (43.1) | 78 (53.8) | 65 (34.8) |
| &nbsp;&nbsp;2 | 43 (13.0) | 26 (17.9) | 17 (9.1) |
| MMSE, mean (SD) | 27.1 (1.8) | 26.7 (1.7) | 27.4 (1.7) |
| ADAS-Cog 13, mean (SD)ᵃ | 18.6 (6.3) | 21.2 (5.3) | 16.5 (6.3) |
| CDR-SB, mean (SD) | 1.6 (0.9) | 1.8 (0.9) | 1.4 (0.8) |
| FAQ, mean (SD)ᵇ | 3.8 (4.4) | 5.4 (4.9) | 2.5 (3.4) |
| CC volume, mm³, mean (SD) | 2618 (453) | 2581 (444) | 2647 (459) |
| Intracranial volume, cm³, mean (SD) | 1585 (168) | 1575 (173) | 1592 (165) |

ᵃ Missing for 3 participants (2 converted, 1 not converted). ᵇ Missing for 2 participants (both not converted). CC volume is the sum of the five FreeSurfer CC subregions.

**Table S1. SNP genotype distribution.**

| SNP | Genotype 0 / 1 / 2, n | Coded allele frequency | HWE exact p |
|---|---|---:|---:|
| rs929708 | 164 / 143 / 25 | 0.291 | 0.505 |
| rs11719939 | 190 / 124 / 18 | 0.241 | 0.766 |
| rs4639533 | 16 / 104 / 212 | 0.795 | 0.500 |
| rs2515029 | 113 / 168 / 51 | 0.407 | 0.426 |

## Coefficient functions

Figure 1 shows the estimated SNP effects on the two shape coordinates, α₁(s) and α₂(s), and the shape coefficient functions in the hazard model, β₁(s) and β₂(s). [TODO: describe where α(s) departs most from zero and in which direction.] REML retained curvature in β₁(s) (effective degrees of freedom [EDF] 3.67–3.73 across the four SNP models) and reduced β₂(s) to a linear function (EDF 2.00). Using 10 or 30 instead of 20 basis functions per coordinate changed the EDF of β₁(s) by at most 0.11 and the indirect HRs by less than 0.001.

**Figure 1.** Estimated coefficient functions α₁(s), α₂(s), β₁(s), and β₂(s) along the CC contour. (`results/penalized_functional_cox/penalized_functional_cox_functions_k20.png`)

## Global tests of the two mediation paths

Conditional on the SNP and covariates, CC shape was associated with the hazard of conversion in all four models (permutation p < 0.001; Table 2). This association came mainly from the first shape coordinate (approximate p < 0.001 for β₁; p = 0.26–0.28 for β₂).

The SNP-to-shape path was weaker. In the integrated test over the whole contour, the association was nominally significant for rs4639533 (p = 0.022) and rs2515029 (p = 0.019), but not for rs929708 (p = 0.114) or rs11719939 (p = 0.212). Neither nominal result remained below 0.05 after Holm correction across the four SNPs (both p = 0.078). Tests based on the maximum pointwise statistic were not significant for any SNP (p ≥ 0.24), indicating that no single contour position carried a strong SNP effect.

**Table 2. Global permutation tests of the two paths (10,000 permutations).**

| SNP | SNP → shape, integrated test p | Holm-adjusted p | SNP → shape, maximum test p | Shape → hazard p |
|---|---:|---:|---:|---:|
| rs929708 | 0.114 | 0.228 | 0.390 | < 0.001 |
| rs11719939 | 0.212 | 0.228 | 0.770 | < 0.001 |
| rs4639533 | 0.022 | 0.078 | 0.240 | < 0.001 |
| rs2515029 | 0.019 | 0.078 | 0.432 | < 0.001 |

SNP → shape: H₀ α₁(s) = α₂(s) = 0 for all s, given the covariates. The integrated test uses the integral of the squared pointwise t statistic over both coordinates, and the maximum test uses the largest absolute pointwise t statistic. Shape → hazard: H₀ β₁(s) = β₂(s) = 0 for all s, given the SNP and the covariates, tested by the partial likelihood ratio of the penalized model, with smoothing parameters reselected in every permutation. Both are Freedman–Lane permutation tests that permute subjects' whole residual curves.

## Direct, indirect, and total effects

Table 3 shows the decomposition of each SNP effect. The indirect effect through CC shape excluded HR = 1 for rs4639533 (HR 1.135, 95% BCa CI 1.018–1.307) and rs929708 (HR 1.095, 95% BCa CI 1.002–1.219). The interval for rs11719939, which is in strong linkage disequilibrium with rs929708, narrowly included 1 (HR 1.091, 95% BCa CI 0.999–1.215). The interval for rs2515029 included 1 (HR 0.950, 95% BCa CI 0.863–1.034). Percentile intervals led to the same conclusions.

rs4639533 had the most consistent evidence for shape mediation: its indirect-effect interval excluded 1, and both of its path tests were below 0.05. For rs929708, the indirect-effect interval excluded 1 by a small margin, but its SNP-to-shape test was not significant.

No direct or total effect excluded HR = 1. For rs4639533 and rs929708, the point estimates of the direct effects (HR 0.864 and 0.840) were opposite in direction to the indirect effects, which partly explains why the total effects were close to 1 (HR 0.981 and 0.920).

Within the indirect effect, the two shape coordinates contributed in opposite directions, and the first coordinate dominated. On the log-HR scale, the contributions of the first and second coordinates were +0.189 and −0.062 for rs4639533, +0.122 and −0.032 for rs929708, +0.113 and −0.026 for rs11719939, and −0.109 and +0.058 for rs2515029.

**Table 3. Direct, indirect, and total effects (multiplier wild bootstrap, 10,000 replicates).**

| SNP | Direct HR (95% BCa CI) | Indirect HR (95% BCa CI) | Indirect, 95% percentile CI | Total HR (95% BCa CI) |
|---|---|---|---|---|
| rs929708 | 0.840 (0.637–1.121) | **1.095 (1.002–1.219)** | 1.001–1.216 | 0.920 (0.687–1.256) |
| rs11719939 | 0.997 (0.753–1.326) | 1.091 (0.999–1.215) | 0.997–1.212 | 1.088 (0.808–1.475) |
| rs4639533 | 0.864 (0.661–1.152) | **1.135 (1.018–1.307)** | 1.014–1.302 | 0.981 (0.729–1.348) |
| rs2515029 | 1.008 (0.793–1.271) | 0.950 (0.863–1.034) | 0.860–1.032 | 0.958 (0.741–1.236) |

HRs are per additional coded allele. Bold: interval excludes 1. The wild bootstrap holds the smoothing parameters at their REML values; see Table 5 for intervals that reselect them.

## Spatial distribution of the indirect effect

Figures 2 and 3 show the spatial indirect effect, SAIE(s) = α₁(s)β₁(s) + α₂(s)β₂(s), along the CC contour. For all four SNPs, the point estimate of SAIE(s) was largest in magnitude near s = 0.08 [TODO: anatomical label of this segment], with signs matching the overall indirect effects. The peak values were 1.65 for rs929708, 1.26 for rs11719939, 2.01 for rs4639533, and −1.35 for rs2515029.

In the point estimates, the segment within ±0.1 of the peak accounted for more than the whole integrated indirect effect, and other parts of the contour partly offset it. However, neither the pointwise 95% intervals nor the 95% simultaneous band excluded zero at any of the 100 contour positions. The overall indirect effect is estimated more precisely than SAIE at any single position, because integration averages estimation error along the contour.

**Figure 2.** SAIE(s) with pointwise 95% intervals and the 95% simultaneous band (pairs bootstrap, 10,000 replicates). (`results/penalized_functional_cox_saie/saie_B10000_seed20260905_k20.png`)

**Figure 3.** SAIE(s) point estimates on the mean CC contour. The circle marks the position of the largest |SAIE(s)|. (`results/penalized_functional_cox_saie/saie_on_mean_contour.png`)

## Shape versus volume

To test whether the shape-mediated effects reflected CC size, total CC volume divided by intracranial volume (ICV) was added as a second, scalar mediator (Table 4). The indirect effects through volume were close to 1 for all SNPs (HR 0.996–1.003), and all of their intervals included 1. The indirect effects through shape were essentially unchanged. For rs4639533, the indirect effect through shape was HR 1.131 (95% BCa CI 1.016–1.301), and the combined indirect effect was HR 1.133 (95% BCa CI 1.015–1.307).

The conclusions were the same when ICV was entered as a covariate instead of dividing volume by it. For rs4639533, the indirect effect through shape was HR 1.134 (95% pairs BCa CI 1.000–1.311). In a model with volume as the only mediator, no joint indirect effect through the five CC subregions excluded 1. Of the 20 subregion-specific indirect effects, only one excluded 1: rs2515029 through the mid-anterior subregion (HR 0.965, 95% BCa CI 0.911–0.997). It would not remain significant after Bonferroni correction for 20 tests.

When the five subregion volumes were entered together with shape, the indirect effects through shape were attenuated, and all of their intervals included 1 (rs4639533: HR 1.103, 95% pairs BCa CI 0.959–1.292). This model adds five correlated mediators and is less precise. Overall, the mediation signal was carried by CC shape rather than CC size.

**Table 4. Indirect effects with shape and total CC volume (/ICV) as joint mediators (multiplier wild bootstrap, 10,000 replicates).**

| SNP | Through shape, HR (95% BCa CI) | Through volume, HR (95% BCa CI) | Combined indirect, HR (95% BCa CI) |
|---|---|---|---|
| rs929708 | **1.092 (1.002–1.219)** | 0.999 (0.971–1.011) | 1.091 (0.996–1.219) |
| rs11719939 | 1.088 (1.000–1.212)ᵃ | 0.999 (0.972–1.013) | 1.087 (0.994–1.213) |
| rs4639533 | **1.131 (1.016–1.301)** | 1.003 (0.990–1.046) | **1.133 (1.015–1.307)** |
| rs2515029 | 0.952 (0.863–1.034) | 0.996 (0.960–1.015) | 0.948 (0.859–1.032) |

ᵃ Lower limit 0.9997, so the interval includes 1. Under the pairs bootstrap, all intervals in this table included 1 (through shape for rs4639533: 0.993–1.317).

## Sensitivity analyses

Table 5 summarizes the sensitivity analyses.

**Smoothing-parameter uncertainty.** The wild bootstrap holds the smoothing parameters fixed. When they were reselected by REML in each replicate of a pairs bootstrap, the standard errors of the indirect effects increased by 11–16%. No indirect-effect interval then excluded 1 (rs4639533: HR 1.135, 95% BCa CI 0.997–1.322; rs929708: HR 1.095, 95% BCa CI 0.985–1.224). The wild-bootstrap basic intervals also included 1 for all SNPs (rs4639533: 0.990–1.271).

**Comparison with the FPCA-based model.** In the functional linear Cox regression model with 18 functional principal components, the indirect effects had the same direction and similar magnitude (rs4639533: HR 1.131, 95% BCa CI 0.966–1.368). The penalized model reduced the standard error of the indirect effect by 16–29% under the pairs bootstrap and by 21–37% under the wild bootstrap.

**Proportional hazards.** Schoenfeld residual tests did not indicate departures from proportional hazards for the SNPs (p = 0.68–0.88) or the fitted shape term (p = 0.66). APOE ε4 violated the assumption (p = 0.005; global test p = 0.07). With a separate baseline hazard for each APOE ε4 genotype, the indirect effects were essentially unchanged.

**Effect modification and basis dimension.** There was no evidence of interaction between the SNP and the fitted shape term (likelihood ratio p = 0.32–0.68). Results were insensitive to the number of basis functions.

**Table 5. Sensitivity analyses of the indirect effect (HR, 95% CI).**

| Analysis | rs929708 | rs11719939 | rs4639533 | rs2515029 |
|---|---|---|---|---|
| Main: wild bootstrap, BCa | **1.095 (1.002–1.219)** | 1.091 (0.999–1.215) | **1.135 (1.018–1.307)** | 0.950 (0.863–1.034) |
| Wild bootstrap, basic | 1.095 (0.986–1.198) | 1.091 (0.982–1.193) | 1.135 (0.990–1.271) | 0.950 (0.875–1.050) |
| Pairs bootstrap with REML reselection, BCa | 1.095 (0.985–1.224) | 1.091 (0.975–1.219) | 1.135 (0.997–1.322) | 0.950 (0.857–1.045) |
| FLCRM with 18 FPCs, pairs BCaᵃ | 1.070 (0.907–1.243) | 1.065 (0.885–1.246) | 1.131 (0.966–1.368) | 0.953 (0.834–1.102) |
| APOE ε4-stratified baseline hazard, point estimate | 1.097 | 1.091 | 1.132 | 0.951 |
| 10 / 30 basis functions, point estimate | 1.095 / 1.095 | 1.091 / 1.091 | 1.136 / 1.135 | 0.950 / 0.950 |

Bold: interval excludes 1. ᵃ 5,000 bootstrap replicates; all other bootstraps used 10,000.

## Summary

CC shape was strongly associated with the hazard of conversion to AD, whereas the SNP-to-shape associations were weak. Under the primary analysis, rs4639533 showed the most consistent evidence of an indirect effect through CC shape, and this effect was not explained by CC volume. Its point estimates were largest near one contour segment, although no individual position was significant. However, the SNP-to-shape association of rs4639533 did not remain significant after correction for the four SNPs tested, and its indirect-effect interval included 1 when smoothing-parameter uncertainty was propagated through the bootstrap.
