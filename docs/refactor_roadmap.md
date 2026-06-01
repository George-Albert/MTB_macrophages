# Refactor Roadmap

## Current State

The repository is in a transition from local exploratory analysis scripts toward a reproducible analysis project. The first hardening pass introduced:

- Shared path utilities in `Codes/_shared/project_setup.R`.
- A package preflight script in `scripts/preflight_packages.R`.
- Workflow order metadata in `config/pipeline_order.yml`.
- Repository-level documentation in `README.md` and `docs/`.

## Main Risks to Address Next

1. Path assumptions remain in most scripts.
   Scripts still use hard-coded `Analyses/Inputs` and `Analyses/Outputs` paths. Each active script should source the shared setup helper and resolve paths from the project root.

2. Active workflow and exploratory code are mixed.
   Files such as `bmrs.R`, `brms_1.R`, `Test_brm_fitting.R`, and `Plot_output_bmrs.R` appear exploratory or model-development oriented. They should either be documented as active workflow steps or moved to an `archive/` or `exploratory/` area in a separate commit.

3. Several scripts write derived artifacts into input folders.
   This is useful for step-to-step handoff, but it blurs raw input versus generated intermediate data. Consider splitting `Inputs/` into raw/external inputs and `Outputs/` or `intermediate/` for generated pipeline handoffs.

4. Script names encode chronology but not always purpose cleanly.
   Keep numeric prefixes for provenance, but eventually standardize names without spaces or typo-like variants.

## Recommended Commit Sequence

1. Repository hygiene and documentation.
   Include `.gitignore`, `README.md`, `docs/`, `config/`, `scripts/preflight_packages.R`, and shared setup utilities.

2. Bootstrap critical scripts.
   Commit path/bootstrap changes for `001_QC_and_dim_reduction.R`, `002_Differential_expression.R`, and `001_Setting_exp_data.R`.

3. Bootstrap the remaining active workflow scripts.
   Refactor only headers/path setup first, without changing analytical logic.

4. Classify exploratory code.
   Move or document experimental Bayesian/eQTL scripts after confirming which are part of the intended final workflow.

5. Extract repeated functions.
   Move reusable plotting, clustering, GO, and Hill-fitting helpers into shared modules only after all scripts run with project-relative paths.

## Proposed Stable Structure

```text
Codes/
  _shared/
    project_setup.R
    plotting_helpers.R
    clustering_helpers.R
    go_helpers.R
  DE_codes/
  EQTL_mapping/
config/
docs/
scripts/
Inputs/     # local-only raw/private inputs
Outputs/    # local-only generated results
```

Longer term, consider renaming `Codes/` to `R/` or `analysis/`, but only after the current GitHub migration is stable.
