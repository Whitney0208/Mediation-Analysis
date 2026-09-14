# Complete ADNI Mediation Analysis Dataset

This directory contains the complete subject-level dataset prepared for the two-layer mediation analysis:

```text
SNP exposure -> corpus callosum shape mediator -> time to AD conversion
```

The data are from the ADNI1 cohort. The ShapeMA files provide the corpus callosum shape, SNPs and subject-level clinical information. The SPLS-Cox files provide the survival time and event indicator. Subjects were linked by ADNI `RID`.

## Source data

| Source | Repository / file | Original content | Use here |
|---|---|---|---|
| ShapeMA | [miyeonyeon/ShapeMA](https://github.com/miyeonyeon/ShapeMA) | 707 subjects, aligned corpus callosum SRVF, four SNPs, clinical information and two PC design matrices | Exposure, mediator and candidate covariates |
| SPLS-Cox | [XuqiaoLi/SPLS-Cox](https://github.com/XuqiaoLi/SPLS-Cox) | 372 MCI subjects after the repository removes row 331 | Survival time, event and raw SPLS clinical fields |

The SPLS-Cox repository removes row 331 because the corresponding shape row is missing. The remaining survival records were matched to ShapeMA subjects by `RID`. The overlap is 333 subjects.

## Complete files

| File | Description |
|---|---|
| `mediation_model_complete.rds` | Complete R list containing all analysis objects and metadata |
| `mediation_model_complete.RData` | Same complete object in `.RData` format |
| `mediation_model_complete_analysis_frame.csv` | 333-row subject-level table with 87 columns |
| `mediation_model_complete_variable_catalog.csv` | Variable roles, sources and notes |
| `mediation_model_matched.rds` | Original matched analysis object used by the first model scripts |

The functional mediator `M` cannot be represented conveniently as ordinary CSV columns, so it is stored in the RDS and RData files. The CSV contains the subject-level fields, SNPs, covariates, survival variables and matching information.

## Object structure

| Object | Dimension | Meaning |
|---|---:|---|
| `subject` | 333 x 49 | Original ShapeMA subject table plus `time`, `event` and matching row identifiers |
| `analysis_frame` | 333 x 87 | Subject-level table combining all preserved fields, SNPs, both PC designs and raw SPLS columns |
| `X` | 333 x 4 | Four SNP exposures coded as additive dosage 0/1/2 |
| `M` | 333 x 100 x 2 | Corpus callosum SRVF shape mediator |
| `W_pc2` | 333 x 8 | Intercept, sex, age, handedness, APOE4, PC1, PC2 and education years |
| `W_pc5` | 333 x 11 | Intercept, sex, age, handedness, APOE4, PC1 through PC5 and education years |
| `spls_clinical_raw` | 333 x 15 | Original SPLS survival and unnamed clinical columns, preserved by position |

The last dimension of `M` has size 2 because the SRVF representation has two shape-coordinate functions. The two coordinates jointly form one multivariate functional mediator. The 100 grid points are positions along the normalized corpus callosum contour, not 100 independent mediators.

## Sample and outcome

| Field | Value |
|---|---:|
| Matched subjects | 333 |
| AD conversion events | 145 |
| Censored observations | 188 |
| Event proportion | 43.5% |
| Survival time range | 154 to 1695 on the source time scale |
| Survival time median | 735 |
| Baseline visit | All `VISCODE = bl` |
| Baseline diagnosis | 332 LMCI and 1 AD |

`event = 1` denotes an observed AD conversion in the SPLS-Cox outcome data. `event = 0` denotes censoring. The survival time is retained on the original SPLS-Cox scale.

The matched table contains one subject, `RID 739`, whose ShapeMA baseline diagnosis is `AD`, while the SPLS-Cox record has `event = 0`. This diagnosis definition must be checked before the final analysis. If the subject was already AD at the survival time origin, the LMCI-to-AD analysis should use the 332 LMCI subjects.

## SNP exposures

The four exposures are stored in `X` and are also available in `analysis_frame` with the prefix `X_`.

| SNP | Dosage 0 | Dosage 1 | Dosage 2 | Missing |
|---|---:|---:|---:|---:|
| `rs929708` | 165 | 143 | 25 | 0 |
| `rs11719939` | 191 | 124 | 18 | 0 |
| `rs4639533` | 16 | 105 | 212 | 0 |
| `rs2515029` | 113 | 168 | 52 | 0 |

Dosage 0, 1 and 2 indicate the number of copies of the coded allele. The source files do not document the allele identity in the prepared object, so results should be described as per one additional coded allele until the allele orientation is verified.

## Candidate covariates

All original ShapeMA subject fields remain in `subject` and `analysis_frame`. The complete object also stores the candidate covariate table and its names in `candidate_covariates` and `candidate_covariate_names`.

The main prepared design is `W_pc2`: `male`, `age_z`, `handedness`, `APOE4`, `PC1`, `PC2` and `education_years`. The alternative design `W_pc5` adds `PC3`, `PC4` and `PC5`. `education_years` is copied from the original `PTEDUCAT` field and is included in both designs. The intercept is included for the first-stage shape model and is not a substantive confounder.

Additional preserved candidate fields include `ICV` (intracranial volume), `SITE`, `PTETHCAT`, `PTRACCAT`, `PTMARRY`, baseline cognitive measures such as `ADAS13`, `MMSE`, `CDRSB` and `FAQ`, and biomarkers such as `FDG`, `ABETA`, `PTAU` and `TAU`.

The additional fields are available for selection but are not automatically confounders. Education is complete in the matched sample and is now part of both prepared adjustment sets. `PTRACCAT` has no variation in the matched data, `PTETHCAT` has very limited variation, and `SITE` has 56 observed site codes. Baseline cognition and biomarkers require a causal decision before adjustment because they may represent disease severity or downstream biology.

## Shape mediator

The mediator is the aligned corpus callosum SRVF object `M`, also called `CC_q` in the source MATLAB file. It has dimensions subject x grid position x shape coordinate:

```r
data <- readRDS("mediation_model_complete.rds")
dim(data$M)
# 333 100 2
```

The two shape coordinates are modeled jointly. In the current functional Cox implementation they are reduced by weighted joint FPCA. The main model retains 18 shape FPCs, explaining 85.96% of the observed shape variation.

For each SNP model, the first-stage coefficient table contains the smoothed and raw coefficient functions for the SNP and every column of the selected `W` matrix. The second-stage Cox coefficient table contains the SNP coefficient, each scalar covariate coefficient, and the coefficients of the retained shape FPC scores. The reconstructed functional Cox coefficient for the shape mediator is stored in the model results and the functional-parameter CSV files.

The ancestry covariates `PC1` through `PC5` in `W_pc2` and `W_pc5` are genetic population-structure variables. They are different from the shape FPCs used to represent `CC_q` in the Cox layer.

## Data quality checks

The matched dataset has no missing values in the four SNPs, `W_pc2`, `W_pc5`, `time`, `event` or `M`. There are no duplicated matched RIDs. The shape array is finite and has dimensions 333 x 100 x 2.

The following checks should be completed before final inference:

- Confirm the baseline diagnosis and survival origin for `RID 739`.
- Verify the coded allele orientation for all four SNPs.
- Decide whether `ICV` belongs to the prespecified adjustment set.
- Decide how a center effect for `SITE` should be handled.
- Keep eligibility (`DX.bl`) separate from ordinary adjustment variables.

## Rebuilding the complete object

Run the following script from the project root:

```bash
Rscript code/05_build_complete_analysis_dataset.R
```

The script reconstructs the complete object from the matched data and the original `clinical.dat`, verifies RID and survival order, and writes the four complete-data outputs in this directory.
