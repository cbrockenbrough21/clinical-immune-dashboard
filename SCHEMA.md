# Database Schema

The schema is normalized to support scalability and analytical flexibility.

--- 

## Modeling Assumptions

- **Treatment assignment is subject-level**: a subject’s treatment is assumed not to vary across their samples.
  - The dataset may include `treatment = 'none'` for subjects without a treatment assignment.
- **Response is subject-level**: response is assumed to be consistent across all samples for a subject.
  - Response may be missing/null for untreated subjects (commonly when `treatment = 'none'`).
  - Ingestion validates that any non-null response values are consistent per subject; conflicting non-null responses fail fast.
- **Project membership**: for this dataset, each subject is assumed to belong to a single project.
  - If subjects can belong to multiple projects in future datasets, the schema can be extended with a `projects` table and a join table (e.g., `subject_projects`).

---

## subjects

Stores subject-level metadata. Each row represents a unique enrolled patient.

| Column | Type | Description |
|--------|------|-------------|
| subject_id | TEXT (PK) | Unique patient identifier |
| project | TEXT | Project identifier |
| condition | TEXT | Clinical condition (e.g., melanoma) |
| age | INTEGER | Age at enrollment |
| sex | TEXT | Biological sex (M/F) |
| treatment | TEXT | Treatment assigned to the subject |
| response | TEXT | Clinical response outcome (yes/no) |

### Modeling Assumptions

- Treatment is modeled at the subject level under the assumption that treatment assignment does not vary across samples for a given subject.
- Response is assumed to be a subject-level outcome and remains consistent across all samples for that subject.
---

## samples

Stores biological sample-level metadata.

| Column | Type | Description |
|--------|------|-------------|
| sample_id | TEXT (PK) | Unique sample identifier |
| subject_id | TEXT (FK) | References subjects.subject_id |
| sample_type | TEXT | e.g., PBMC |
| time_from_treatment_start | INTEGER | Timepoint (days) |

---

## cell_counts

Stores immune population measurements in long format.

| Column | Type | Description |
|--------|------|-------------|
| sample_id | TEXT (FK) | References samples.sample_id |
| population | TEXT | Immune cell type |
| count | INTEGER | Cell count |

Primary Key: (sample_id, population)

---

## Indexes

- idx_subjects_project ON subjects(project)
- idx_subjects_treatment ON subjects(treatment)
- idx_samples_subject ON samples(subject_id)
- idx_samples_type_time ON samples(sample_type, time_from_treatment_start)
- idx_cell_counts_population ON cell_counts(population)

---

## Design Rationale

- Separates subject-level and sample-level concepts.
- Uses long format for immune populations to allow extensibility.
- Supports indexing for fast filtering by condition, treatment, or response.
- Enables efficient aggregation queries for statistical analysis.