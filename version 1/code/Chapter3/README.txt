Chapter 3 local reproduction files

Main script:
- chapter3_reproduce_local.R

Local function copies used by the script:
- fn/1.R
- fn/3.R
- fn/5.R
- fn/7.R

Note:
- chapter3_reproduce_local.R now sources local helper functions from `Chapter3/fn`
- it does not rely on editing the original project files in `sgmas - anisha - project3`

What it does:
- reads Q3_Covariates.csv and LQD_transformed_pixel.RDS from:
  /Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/input_data/Chpater3
- runs the survHIMA branch
- runs the hima2 survival branch
- writes outputs into:
  /Users/wangwanbing/Desktop/Chao Lab/Mediation Analysis/wwanbing/output/Chapter3

Outputs:
- chapter3_mediation_survHIMA.csv
- chapter3_mediation_hima2.csv
- chapter3_mediation_survHIMA_ide.png
- chapter3_mediation_hima2_ide.png
