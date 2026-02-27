# Clinical Immune Dashboard

A reproducible clinical immune profiling pipeline built with Python and SQLite.

This project ingests immune cell population data, computes sample-level frequency metrics, performs responder stratification analysis, and exposes results via a lightweight dashboard.

---

## Overview

This system implements:

### Data Ingestion Layer
- Normalized relational schema (SQLite)
- Deterministic database rebuilds
- Validation of subject-level response consistency
- Long-format immune population storage for scalability

### Sample-Level Immune Frequency Metrics
- Total cell counts per sample
- Relative immune population frequencies (%)
- Exported analytical dataset for downstream analysis

### Responder Stratification Analysis
- Melanoma PBMC cohort
- Miraclib-treated patients
- Mann–Whitney U test per immune population
- Benjamini–Hochberg false discovery rate correction
- Publication-ready boxplots

### Baseline Cohort Exploration
- Time-zero melanoma PBMC subset
- Cohort summaries by project, response, and sex
- Baseline B-cell summary statistic for melanoma male responders

---

## Repository Structure

```
clinical-immune-dashboard/
│
├── load_data.py                     # Database initialization and ingestion
├── run_analysis.py                  # Analytical workflows
├── Makefile                         # setup / pipeline / dashboard targets
├── requirements.txt
│
├── data/
│   └── cell-count.csv               # Source dataset
│
├── outputs/                         # Generated analytical artifacts
│   ├── sample_population_frequencies.csv
│   ├── responder_differential_analysis.csv
│   ├── responder_stratification_boxplot.png
│   ├── baseline_cohort_summary.csv
│   └── baseline_bcell_summary.txt
│
├── SCHEMA.md                        # Database schema and modeling rationale
└── README.md
```

---

## Running the Project

### Install dependencies
```
make setup
```

### Execute full analytical pipeline
```
make pipeline
```

This command will:

1. Initialize the SQLite database (`immune_trial.db`)
2. Load all rows from `data/cell-count.csv`
3. Compute sample-level immune population frequencies
4. Perform responder stratification analysis
5. Generate plots and baseline cohort summaries
6. Write all outputs to the `outputs/` directory

### Launch the dashboard
```
make dashboard
```

The dashboard runs locally (URL will be displayed in the terminal once implemented).

---

## Schema and Modeling Notes

The relational schema is documented in `SCHEMA.md`. The design separates:

- Subjects (patient-level metadata)
- Samples (biological sample-level metadata)
- Cell counts (immune populations in long format)

### Modeling Assumptions

- **Input location**: the pipeline reads from `data/cell-count.csv` (committed to this repo for reproducibility).
- **Treatment is subject-level**: treatment assignment is assumed not to vary across samples for a given subject.
  - The dataset may include `treatment = 'none'` to indicate no treatment assignment.
- **Response is subject-level**: response is assumed to be consistent across all samples for a subject.
  - Response may be missing/null for untreated subjects (e.g., `treatment = 'none'`).
  - During ingestion, we validate that any *non-null* response values for the same subject are consistent.
- **Project membership**: for this dataset, each subject is assumed to belong to a single project.
  - The schema is written to be extensible to multiple projects (e.g., adding a `projects` table).
- **Immune populations are stored in long format** to allow extension to additional cell types without schema changes.
---

## Statistical Methods

Responder vs non-responder comparisons use:

- Mann–Whitney U test (non-parametric)
- Benjamini–Hochberg correction for multiple hypothesis testing (FDR control)

These methods are appropriate for non-normally distributed immune cell frequency data.

---

## Scalability Considerations

The schema and analytical design support:

- Hundreds of projects
- Thousands of subjects
- Millions of immune population measurements
- Additional downstream analytics without schema refactoring

Indexes are defined on primary join keys and common filtering columns to ensure efficient cohort selection as data scales.

---

## Reproducibility

The project is designed to be reproducible from scratch:

- `make setup` installs all dependencies
- `make pipeline` rebuilds the database and regenerates all outputs
- No manual steps are required
