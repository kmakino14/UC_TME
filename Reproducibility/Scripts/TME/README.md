# TME bulk projection

This repository reproduces TCGA bulk RNA-seq TME projection analyses from archived Scaden prediction matrices.

## Data availability

The archived input data for TME deconvolution are not stored directly in this Git repository.

Please download `TME_deconvolution.zip` from the GitHub Release:

`v1.0-tme-deconvolution`

The archive contains the following two directories at its top level:

```text
bulk/
predictions/
```

After downloading the archive, extract it into data/archived_inputs/:

```bash
mkdir -p data/archived_inputs
unzip TME_deconvolution.zip -d data/archived_inputs
```

After extraction, the repository should contain:

```text
data/archived_inputs/bulk/
data/archived_inputs/predictions/
```

The trained Scaden M256 model used to generate these matrices is included for transparency:

```text
model/m256/
```

Model-based re-prediction is not part of the main workflow. The downstream analyses are reproduced directly from the archived prediction matrices.

## Preparing the archived input data

Run from the repository root:

```bash
unzip TME_deconvolution.zip
```

This should create the `data/` directory required by the workflow. If the archive extracts into a nested directory, move the extracted `data/` directory to the repository root before running the scripts.

The expected directory structure is:

```text
data/
└── archived_inputs/
    ├── bulk/
    └── predictions/
```

## Workflow

Run from the repository root after preparing the `data/` directory:

```bash
bash scripts/00_use_archived_inputs.sh
bash scripts/01_run_TME_projection.sh
bash scripts/02_run_TCGA_Cox.sh
bash scripts/03_run_assignment_heatmap.sh
```
