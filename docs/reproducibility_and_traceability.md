# Reproducibility and Traceability

## Reproducibility Controls

- Shared bootstrap utility at `Codes/_shared/project_setup.R`.
- Path resolution supports both layouts:
  - `Inputs/` and `Outputs/`
  - `Analyses/Inputs/` and `Analyses/Outputs/` (legacy compatibility)
- Package preflight script at `scripts/preflight_packages.R`.
- Parse validation script at `scripts/validate_r_parse.R`.
- Pipeline execution order documented in `config/pipeline_order.yml`.

## Traceability Controls

- Numeric script IDs (`001_`, `002_`, etc.) represent execution chronology.
- Inputs and outputs are path-derived from script IDs to maintain provenance.
- Legacy/trash content is excluded in `.gitignore` to prevent accidental publication.

## Recommended Next Hardening Steps

1. Add script-level parameter files per workflow step (`config/*.yml`).
2. Externalize thresholds currently embedded in code.
3. Add deterministic seeds where stochastic routines are used.
4. Add automated smoke tests for key scripts using reduced test inputs.
