# mels-bayes — Bayesian mixed-effects location-scale EMA analysis

Analysis code for: *Behavioral predictors of mood location and scale in depressed youth:
A Bayesian mixed-effects location-scale EMA analysis* (submitted to BMC Medicine).

## What this repository contains

- `scripts/modules/`, `scripts/pipeline/`  模型拟合与主管线（入口: `scripts/pipeline/master.R`）
- `scripts/sensitivity/`  敏感性分析（边界删失、先验、horseshoe 等）
- `scripts/tools/`  图件重导与诊断工具
- `scripts/validation/`  验证脚本

## Reproducing

1. Install R (see `reproducibility.md` for versions) with brms/cmdstanr (CmdStan backend).
2. Run scripts from the repository root: `Rscript scripts/pipeline/master.R`.
3. The de-identified dataset is published separately (see Data below);
   raw data are not redistributed in this repository.

## Data

The de-identified EMA and pre/post dataset (N = 236) with codebook and
de-identification README is deposited on ScienceDB (DOI: to be filled upon publication).
The raw-data cleaning pipeline is not part of this repository: it processes raw,
identified data (never redistributed) and cannot run without it; the de-identified
dataset above is the analysis-ready form.

## License & citation

Code: MIT (see `LICENSE`). If you reuse the code or the dataset, please cite the
paper and the dataset DOI.

