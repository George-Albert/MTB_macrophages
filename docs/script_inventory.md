# Script Inventory

This inventory summarizes the current role of the tracked scripts and the main
technical issues observed during repository hardening.

## Differential Expression Workflow

| Script | Current role | Review status |
| --- | --- | --- |
| `Codes/DE_codes/001_QC_and_dim_reduction.R` | QC, sample/metadata alignment, PCA/t-SNE, processed expression outputs for downstream analyses. | Bootstrap started; still long and monolithic. |
| `Codes/DE_codes/002_Differential_expression.R` | `edgeR`/`limma` differential expression across infection and time contrasts. | Bootstrap started; analytical logic preserved. |
| `Codes/DE_codes/003_Clustering.R` | Selects infection-responsive genes, normalizes temporal logFC profiles, estimates k-means/hierarchical k-means clusters, writes clustering intermediates. | Needs shared bootstrap; contains reusable plotting/clustering helpers that could move to `_shared`. |
| `Codes/DE_codes/004_Clustering_2.R` | Downstream inspection and visualization of selected cluster assignments. | Needs shared bootstrap. |
| `Codes/DE_codes/005_1Setting_data_annotation.R` | Builds GO annotation lists from ClueGO/GO tables and maps identifiers. | Needs shared bootstrap; writes generated annotation intermediates. |
| `Codes/DE_codes/005_Boolean_matrix_building.R` | Builds GO boolean matrices. | Needs shared bootstrap. |
| `Codes/DE_codes/006_Metadata_matrix_builder.R` | Builds GO metadata matrices. | Needs shared bootstrap. |
| `Codes/DE_codes/007_Significant_responses_over_ time.R` | Summarizes significant response patterns over time. | Needs shared bootstrap; filename contains an internal space. |
| `Codes/DE_codes/008_Enrichment_Analysis_by_cluster.R` | Performs GO enrichment by cluster/trend and writes enrichment tables. | Needs shared bootstrap; several reusable functions could be extracted later. |
| `Codes/DE_codes/009_Enrichment_Heatmaps_by_stage.R` | Builds enrichment heatmaps by response stage. | Needs shared bootstrap. |
| `Codes/DE_codes/010_cor_pvalues_GO.R` | Correlation/p-value analysis over GO outputs. | Needs shared bootstrap. |

## Individual Response and eQTL-Oriented Workflow

| Script | Current role | Review status |
| --- | --- | --- |
| `Codes/EQTL_mapping/001_Setting_exp_data.R` | Builds individual-level TB-minus-NI response matrices across time and joins cluster labels. | Bootstrap started; one remaining clustering path still uses legacy `Analyses/Inputs`. |
| `Codes/EQTL_mapping/002_Fitting_Hill_eq.R` | Fits nonlinear Hill-type curves to temporal response profiles by gene, individual, and cluster. | Needs shared bootstrap; currently appears configured to run only the gradual-up cluster set. |
| `Codes/EQTL_mapping/003_Plot_Hill_fits_wo_constrains.R` | Plots Hill fits without constraints. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/004_Inspect_R2_negative.R` | Inspects fits with negative R-squared values. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/005_Test_filtering_options_2.R` | Tests filtering thresholds over Hill-fit quality metrics and individual coverage. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/006_Random_inspection_w_constrains.R` | Inspects constrained fits and parameter distributions. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/008_Correlation_parameters.R` | Correlates Hill model parameters and prepares selected parameter groups. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/009_R2_min_Ind.R` | Evaluates R2 and minimum-individual thresholds. | Needs shared bootstrap. |
| `Codes/EQTL_mapping/bmrs.R`, `brms_1.R`, `Test_brm_fitting.R`, `Plot_output_bmrs.R` | Bayesian/nonlinear model exploration and plotting. | Should be classified as active workflow or moved to exploratory/archive in a separate commit. |
| `Codes/EQTL_mapping/Raul_code/run.01.matrixeQTL.R` | External/legacy Matrix eQTL-oriented workflow. | Keep documented separately unless it becomes part of the active pipeline. |

## Cross-Cutting Technical Notes

- Many scripts still assume `Analyses/Inputs` and `Analyses/Outputs`.
- Several scripts write generated intermediates into input directories; this
  should be documented or gradually moved to an explicit intermediate area.
- The main analytical scripts are long and mix setup, modelling, plotting, and
  output writing. Function extraction should happen only after path migration.
- The first safe refactor pass should change script headers/path setup only,
  without changing statistical models or thresholds.
