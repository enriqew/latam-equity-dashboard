# latam-equity-dashboard

Pipeline and data engineering source for the **LATAM Equity Dashboard** — a visualization of the double inequity faced by Latin American patients in organ transplantation.

## Concept: Double Inequity

The dashboard synthesizes two independent dimensions of transplant inequity:

| Dimension | Metric | Source |
|---|---|---|
| **Access inequity** | Transplants per million population (pmp) | IRODaT / GODT |
| **Pharmacogenomic inequity** | Delta vs CEU baseline (pp) for CPIC A-level drugs | 1000 Genomes × CPIC |

A patient or population in the **top-left quadrant** of the scatter plot faces both fewer transplant opportunities *and* a higher probability that standard European-calibrated immunosuppression protocols will be sub-optimal for their genotype.

## Data Inputs

### transplant-atlas (access dimension)

| File | Fields | Description |
|---|---|---|
| `world-transplants.json` | `country_iso3`, `country_name`, `year`, `total_transplants`, `transplants_pmp`, `deceased_donors_pmp`, `living_donors_pmp`, `source` | IRODaT global transplant volumes, 104 countries, 1993–2025 |

### pgx-latam-atlas (pharmacogenomic dimension)

| File | Fields | Description |
|---|---|---|
| `drug_impact_summary.json` | `drug_name`, `gene_symbol`, `population_code`, `percentage_requiring_change`, `baseline_ceu_percentage`, `delta_vs_baseline` | % of individuals per population requiring dose/therapy change per CPIC guidelines |
| `actionability_ranking.json` | `drug_name`, `gene_symbol`, `population_code`, `delta_vs_baseline`, `population_affected_pct` | Top-ranked drug–gene–population combinations by clinical impact vs CEU baseline |

## Population → Country Mapping

The join between PGx and transplant data uses a hardcoded mapping from 1000 Genomes superpopulation codes to country ISO3:

| Population Code | Country ISO3 | Country |
|---|---|---|
| MXL | MEX | Mexico |
| PEL | PER | Peru |
| CLM | COL | Colombia |
| PUR | PRI | Puerto Rico |
| CEU | GBR | United Kingdom (European proxy) |

### Limitations of this mapping

1. **Admixed proxies, not national surveys.** MXL (Mexican Ancestry from LA, CA), PEL (Peruvians from Lima), CLM (Colombians from Medellin), and PUR (Puerto Ricans from Puerto Rico) are 1000 Genomes samples — they represent urban, admixed subpopulations, not exhaustive national allele frequency surveys.

2. **Puerto Rico / IRODaT gap.** IRODaT does not report Puerto Rico as a separate jurisdiction. The dashboard uses an estimate derived from US territory data (~45 pmp). This is flagged explicitly in the UI.

3. **CEU → GBR proxy.** CEU (Utah residents with Northern and Western European ancestry) is used as the European baseline. GBR transplant data is used as the access comparator for this population.

4. **Colombia 2025 data.** IRODaT 2025 Colombia records lack `transplants_pmp` (total volume only). The dashboard falls back to the most recent year with a computed pmp value.

## Equity Score Formula

The double inequity score is conceptual and not computed as a single scalar in the current version. The scatter plot encodes it visually:

```
x-axis: transplants_pmp        (higher = better access)
y-axis: delta_vs_baseline (pp) (higher = larger PGx divergence from European protocols)

Compounded inequity = low x AND high y (top-left quadrant)
```

A future version may compute a composite index:

```python
equity_score = (1 / transplants_pmp_normalized) * (1 + pgx_delta_normalized)
# where both dimensions are normalized 0–1 within the dataset
```

## Planned Pipeline Stack

| Layer | Tool | Purpose |
|---|---|---|
| Ingestion | Python + boto3 | Fetch IRODaT data, 1000 Genomes VCFs |
| Storage | AWS S3 + Apache Iceberg | Raw and processed data lake |
| Transformation | AWS Glue + dbt | Clean, normalize, compute deltas |
| Query | Athena + DuckDB | Ad-hoc analysis and validation |
| Export | Python (pandas → JSON) | Static JSON artefacts for portfolio |

## Output Artefacts

The pipeline produces static JSON files committed to the portfolio repo (`data-dive-design-hub/src/data/`):

- `src/data/transplant-atlas/world-transplants.json` — from transplant-atlas pipeline
- `src/data/pgx/drug_impact_summary.json` — from pgx-latam-atlas pipeline
- `src/data/pgx/actionability_ranking.json` — from pgx-latam-atlas pipeline

No runtime API calls. The React dashboard imports these files directly at build time.

## References

- IRODaT: International Registry in Organ Donation and Transplantation — https://www.irodat.org
- CPIC: Clinical Pharmacogenomics Implementation Consortium — https://cpicpgx.org
- 1000 Genomes Project Phase 3 — https://www.internationalgenome.org
