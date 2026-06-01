# MTB Human Macrophages Dynamics

Analysis repository for time-resolved host response to *Mycobacterium tuberculosis*
infection in human macrophages.

## Scientific Context

This project supports a doctoral research line focused on developing
bio-computational methods to study host-pathogen interactions between
*Mycobacterium tuberculosis* and macrophages.

The core biological setting is an in vitro infection experiment in which
macrophages derived from approximately 100 individuals were infected with a
virulent *M. tuberculosis* strain. Host and bacterial RNA were profiled across a
time course after infection, with non-infected controls kept in parallel.

The broader research objective is to move beyond static eQTL analyses and model
how genetic variation influences dynamic transcriptional responses to infection.
The current repository focuses on the host transcriptomic workflow: quality
control, differential expression, temporal response clustering, functional
enrichment, and early modelling of individual response dynamics with nonlinear
Hill-type curves.

## Current Analysis Scope

The active code currently covers:

- RNA-seq quality control and dimensionality reduction.
- Differential expression of infected versus non-infected macrophages across
  time points.
- Selection and clustering of infection-responsive genes based on temporal
  log-fold-change profiles.
- GO/functional enrichment and enrichment heatmap generation.
- Construction of individual-level infection response matrices.
- Hill-equation fitting of temporal response profiles for downstream
  QTL-oriented analyses.

The repository is being professionalized while preserving analytical behaviour.
Large private inputs and generated outputs remain local-only.

## Repository Structure

- `Codes/DE_codes/`: Differential expression, clustering, and enrichment workflow.
- `Codes/EQTL_mapping/`: Individual response modelling and eQTL-oriented scripts.
- `Codes/_shared/`: Shared path/bootstrap utilities for reproducible execution.
- `config/`: Pipeline metadata and execution order.
- `scripts/`: Utility scripts for preflight and reproducibility checks.
- `docs/`: Traceability, reproducibility, review notes, and refactor roadmap.
- `Inputs/`, `Outputs/`: Local data/artifacts, intentionally excluded from git.

## Workflow Overview

The recommended execution order is documented in `config/pipeline_order.yml`.

At a high level:

1. `001_QC_and_dim_reduction.R` prepares expression and metadata tables after QC.
2. `002_Differential_expression.R` models host expression with `edgeR`/`limma`.
3. `003_Clustering.R` clusters infection-responsive genes by temporal response.
4. `004_Clustering_2.R` refines and visualizes selected cluster structure.
5. `005_*` to `010_*` build GO resources, enrichment summaries, and heatmaps.
6. `EQTL_mapping/001_Setting_exp_data.R` builds individual infection response
   matrices and joins cluster labels.
7. `EQTL_mapping/002_*` onward fit and inspect nonlinear temporal response models.

See `docs/script_inventory.md` for a current review of script roles and
remaining hardening needs.

## Reproducible Execution

Scripts are being migrated to resolve paths through the shared bootstrap helper:

```r
source("Codes/_shared/project_setup.R")
paths <- bootstrap_project()
```

This resolves:

- the project root;
- `Inputs/` or legacy `Analyses/Inputs/`;
- `Outputs/` or legacy `Analyses/Outputs/`.

Run package preflight before executing analysis scripts:

```r
source("scripts/preflight_packages.R")
```

On this machine, `Rscript` was found at:

```text
C:/Users/JorgeAlbertoCardenas/AppData/Local/Programs/R/R-4.6.0/bin/Rscript.exe
```

## Data Policy

This repository tracks code, configuration, and documentation. It does not track:

- raw RNA-seq/genotype inputs;
- processed expression matrices;
- generated RDS/RData/model outputs;
- local planning PDFs;
- trash/archive folders and editor/session artifacts.

These exclusions are intentional because the local data include large files and
potentially private research artifacts.

## Refactor Status

The repository is being hardened in conservative passes:

- preserve analytical logic;
- centralize path and package handling;
- document workflow order and generated artifacts;
- separate active workflow scripts from exploratory/model-development scripts;
- gradually migrate remaining scripts away from hard-coded `Analyses/...` paths.

See `docs/refactor_roadmap.md` for the current restructuring plan.
