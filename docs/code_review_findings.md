# Code Review Findings

## High Priority

1. Path assumptions still present in many scripts.
   - Most scripts still expect specific relative directories and rely on the current working directory.
   - Mitigation started by introducing `Codes/_shared/project_setup.R` and updating three critical scripts.
   - Remaining high-impact candidates: `003_Clustering.R`, `004_Clustering_2.R`, `005_1Setting_data_annotation.R`, `008_Enrichment_Analysis_by_cluster.R`, and EQTL scripts `002_` through `009_`.

2. Monolithic scripts with mixed responsibilities.
   - Several files combine data preparation, analysis, and plotting in one script.
   - Risk: harder debugging, lower reusability, and lower testability.

## Medium Priority

1. Package loading is duplicated across scripts.
   - Increases drift risk when dependencies change.
   - Preflight script centralizes validation but load blocks are still repeated.

2. Inconsistent comment language and style.
   - Existing comments mix Spanish and English and vary in granularity.
   - Initial translation to English added in the first DE script section.

## Low Priority

1. Legacy/experimental scripts coexist with active workflow scripts.
   - Some are already ignored, but naming conventions could be clearer.

2. Minor naming issues.
   - Example: `007_Significant_responses_over_ time.R` contains an internal whitespace artifact.

## Changes Applied In This Revision

- Added shared setup utility: `Codes/_shared/project_setup.R`.
- Added `bootstrap_project()` helper for future script header simplification.
- Refactored headers/path bootstrap for:
  - `Codes/DE_codes/001_QC_and_dim_reduction.R`
  - `Codes/DE_codes/002_Differential_expression.R`
  - `Codes/DE_codes/003_Clustering.R`
  - `Codes/EQTL_mapping/001_Setting_exp_data.R`
  - `Codes/EQTL_mapping/002_Fitting_Hill_eq.R`
- Added professional repo docs and workflow metadata:
  - `README.md`
  - `config/pipeline_order.yml`
  - `docs/reproducibility_and_traceability.md`
- Added refactor roadmap:
  - `docs/refactor_roadmap.md`
- Added script inventory:
  - `docs/script_inventory.md`
- Expanded README scientific context from the local doctoral research plan while
  keeping the planning PDF excluded from git.
