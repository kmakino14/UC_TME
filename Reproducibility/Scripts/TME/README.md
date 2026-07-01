# UC multiome TME bulk projection

This repository reproduces TCGA bulk RNA-seq TME projection analyses from archived Scaden prediction matrices.

## Overview

The workflow starts from archived bulk input files and archived Scaden prediction matrices:

```text
data/archived_inputs/bulk/
data/archived_inputs/predictions/
```

The trained Scaden M256 model used to generate these matrices is included for transparency:

```text
model/m256/
```

Model-based re-prediction is not part of the main workflow. The downstream analyses are reproduced directly from the archived prediction matrices.

## Workflow

Run from the repository root:

```bash
bash scripts/00_use_archived_inputs.sh

conda activate R4_3
bash scripts/01_run_TME_projection.sh
bash scripts/02_run_TCGA_Cox.sh
bash scripts/03_run_assignment_heatmap.sh
```

## Main scripts

```text
scripts/00_use_archived_inputs.sh
scripts/01_run_TME_projection.sh
scripts/02_run_TCGA_Cox.sh
scripts/03_run_assignment_heatmap.sh
```

## Core analysis scripts

```text
scripts/core/TME_projection.R
scripts/core/run_TME_projection.sh
scripts/core/TCGA_multivariable_Cox.R
scripts/core/assignment_heatmap.R
scripts/core/run_assignment_heatmap.sh
```

## Default output

```text
/home/k-makino/wd_home/UC_DOGMA_reseq/Reproducibility/Results/TME
```

The output structure is:

```text
TME/
├── BLCA/archived_m256/
│   ├── bulk/
│   ├── pred/
│   └── module_signature_template_projection/
├── BRCA/archived_m256/
├── ...
├── log/
├── TCGA_multivariable_Cox/
└── TCGA_assignment_evidence_heatmap/
```

## TME categories

Final bulk TME assignment uses:

```text
TME1, TME3, TME4, TME5, TME6
```

The excluded category is documented here for interpretation only and is not included in output directory or file names.

## Requirements

The downstream R scripts were tested in the `R4_3` conda environment.
