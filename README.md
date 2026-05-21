# latam-equity-dashboard

A Python data pipeline that computes a **dual-inequity index** for kidney transplant patients in Latin America. The index combines two independent dimensions of disadvantage: access to transplantation (how many transplants happen relative to population size) and pharmacogenomic divergence (how much LATAM populations differ from the European baseline that drives most clinical dosing guidelines).

The output is a set of JSON artifacts ready for visualization or downstream analysis. There is no machine learning and no database; the pipeline is a functional composition of ingest, transform, and export steps.

---

## Table of Contents

1. [Background and Motivation](#background-and-motivation)
2. [Index Methodology](#index-methodology)
   - [Original Index — 4 Countries](#original-index--4-countries)
   - [Extended Index — 18 Countries](#extended-index--18-countries)
   - [Score Formula](#score-formula)
3. [CPIC Gene-Drug Pairs](#cpic-gene-drug-pairs)
4. [Exclusions and Rationale](#exclusions-and-rationale)
5. [Data Sources](#data-sources)
6. [Puerto Rico Note](#puerto-rico-note)
7. [gnomAD Hardy-Weinberg Methodology](#gnomad-hardy-weinberg-methodology)
8. [Project Structure](#project-structure)
9. [Installation](#installation)
10. [Running the Pipeline](#running-the-pipeline)
11. [Artifact Schema](#artifact-schema)
12. [Running Tests](#running-tests)
13. [Clinical Disclaimer](#clinical-disclaimer)

---

## Background and Motivation

Kidney transplantation is the most effective treatment for end-stage renal disease, yet access rates across Latin America vary by more than an order of magnitude relative to high-income countries. At the same time, most pharmacogenomic dosing guidelines — including those published by the Clinical Pharmacogenomics Implementation Consortium (CPIC) — were calibrated on European-ancestry cohorts. Patients from admixed Latin American populations carry different allele frequencies for key drug-metabolizing genes, meaning standard dosing protocols may be systematically miscalibrated for them.

This pipeline makes both gaps quantitative and comparable on a single numeric scale, allowing policy-level comparison across countries and over time.

---

## Index Methodology

### Original Index — 4 Countries

Covers Mexico, Peru, Colombia, and Puerto Rico. Pharmacogenomic data come from the 1000 Genomes Project phase 3 populations with a direct geographic mapping:

| 1000G Population | Country | ISO3 |
|-----------------|---------|------|
| MXL | Mexico | MEX |
| PEL | Peru | PER |
| CLM | Colombia | COL |
| PUR | Puerto Rico | PRI |

The PGx gap for each country is the **mean of absolute deltas** across five CPIC Level A gene-drug pairs (see [CPIC Gene-Drug Pairs](#cpic-gene-drug-pairs)):

```
pgx_gap = mean( |delta_vs_baseline| for each of the 5 focal pairs )
```

where `delta_vs_baseline` is the difference in the percentage of patients requiring a dose or therapy change between the LATAM population and the European CEU baseline.

### Extended Index — 18 Countries

Expands coverage to 18 countries by adding 14 additional countries that lack a direct 1000 Genomes population match. These countries use the gnomAD v3 AMR (admixed American) aggregate as a pharmacogenomic proxy.

The extended index uses a single, methodologically consistent metric for all 18 countries: CYP3A5\*3 expressor probability derived from allele frequencies via Hardy-Weinberg equilibrium.

**1000G-mapped countries (4):** MEX, PER, COL, PRI — same geographic mapping as the original index, re-computed under HWE.

**gnomAD AMR-proxy countries (14):**

| Subregion | Countries |
|-----------|-----------|
| Southern Cone | ARG, BRA, CHL, URY |
| Andean / Gran Chaco | BOL, ECU, PRY, VEN |
| Central America | CRI, GTM, HND, NIC, PAN, SLV |

Countries using the AMR proxy share identical PGx values; they are differentiated only by their transplant access figures. Each AMR-proxy record is flagged with `pgx_is_proxy: true` and `pgx_data_source: "gnomad_r3_hwe_amr_proxy"`.

### Score Formula

Both indexes use the same additive formula:

```
access_gap   = global_median_pmp - country_pmp   # positive → below median (worse access)
pgx_gap      = |delta_vs_ceu|                     # always non-negative
equity_score = access_gap + pgx_gap               # higher → worse compounded inequity
```

`global_median_pmp` is the median of transplants-per-million-population across all countries in the IRODaT/GODT world transplants dataset for the most recent available year. It is recomputed at runtime from the ingested data, not hardcoded.

Countries **above** the global median have a negative `access_gap`, which pulls their equity score down relative to countries below the median. The index is intentionally directional: a high score signals compounded disadvantage on both dimensions simultaneously.

---

## CPIC Gene-Drug Pairs

The five focal pairs used in the original composite gap are all CPIC Level A (highest evidence tier), selected for clinical relevance in the kidney transplant setting:

| Gene | Drug | Clinical Role |
|------|------|---------------|
| CYP3A5 | tacrolimus | Primary calcineurin-inhibitor immunosuppressant post-transplant |
| CYP2C19 | clopidogrel | Antiplatelet agent for post-transplant cardiovascular events |
| VKORC1 | warfarin | Anticoagulation management post-transplant |
| DPYD | fluorouracil | Chemotherapy for malignancies in immunosuppressed patients |
| CYP2C9 | fluvastatin | Statin for cardiovascular protection post-transplant |

Each pair contributes one `|delta_vs_baseline|` value to the composite. The composite is the unweighted arithmetic mean across all pairs present in the data.

---

## Exclusions and Rationale

Two genes were evaluated and explicitly excluded from the composite:

**TPMT** — `delta_vs_baseline = 0` for all four LATAM populations. TPMT allele frequencies in MXL, PEL, CLM, and PUR are indistinguishable from the CEU baseline within the 1000 Genomes dataset. Including it would add noise without informative signal and would artificially dilute the composite mean.

**G6PD** — The CEU baseline in the source CPIC drug impact data appears as approximately 100% actionable, while LATAM populations show 1–8%. This near-zero CEU reference value is inconsistent with the known epidemiology of G6PD deficiency (which is rare in European populations) and indicates a pipeline artifact in the upstream data source. Including G6PD would cause it to dominate the composite and distort the index in a way that does not reflect real population pharmacogenomics.

Both exclusions are documented in `src/latam_equity/pgx_composite.py`.

---

## Data Sources

| Source | What it provides | How it is ingested |
|--------|-----------------|-------------------|
| **IRODaT / GODT** (Global Observatory on Donation and Transplantation) | Transplants per million population by country and year | `ingest/transplant_artifacts.py` — copies from a local sibling repo or fetches from a GitHub release of `enriqew/transplant-atlas` |
| **CPIC drug impact summary** | Per-population percentages requiring dose/therapy change, `delta_vs_baseline` for each gene-drug pair | `ingest/pgx_artifacts.py` — copies from a local sibling repo or fetches from a GitHub release of `enriqew/pgx-latam-atlas` |
| **gnomAD v3** (Broad Institute) | CYP3A5\*3 allele counts and allele numbers for 1000G sub-populations and the AMR aggregate | `ingest/gnomad_pgx.py` — live GraphQL query to the gnomAD API at `https://gnomad.broadinstitute.org/api` |

### Ingest source options

The PGx and transplant ingest modules support two source modes controlled by the `--source` flag:

- `local` (default): reads from a local checkout of the upstream pipeline repos at a configured path
- `github`: fetches JSON files from the `main` branch of the respective GitHub repositories

```bash
python -m ingest.pgx_artifacts --source github
python -m ingest.transplant_artifacts --source github
```

The gnomAD ingest always queries the live API.

---

## Puerto Rico Note

IRODaT does not report Puerto Rico as a separate jurisdiction; transplant volume is aggregated into United States totals. When no Puerto Rico row is found in the transplant data, the pipeline falls back to a fixed estimate of **45 transplants per million population**, derived from US territory-level data. This estimate is applied only to PRI and is flagged in every output record:

```json
"transplants_pmp": 45.0,
"transplants_pmp_is_estimate": true
```

Consumers of the artifact should treat PRI access figures with corresponding caution.

---

## gnomAD Hardy-Weinberg Methodology

For the extended index, CYP3A5 metabolizer status is inferred from allele frequency data using Hardy-Weinberg equilibrium.

**Variant:** CYP3A5\*3 (rs776746), GRCh38 position `7-99672916-T-C`

- The ALT allele (C on the forward strand) corresponds to the \*3 non-functional allele.
- Individuals carrying at least one \*1 (REF T) allele are considered **expressors**: they have functional CYP3A5 and higher tacrolimus clearance, typically requiring higher doses.

**Calculation:**

```
AF(*3)          = ac / an                         # allele frequency from gnomAD
p_expressor     = 1 - AF(*3)^2                    # P(≥1 functional *1 allele) under HWE
pct_expressor   = p_expressor × 100               # percentage of population

delta_vs_ceu_pp = (p_expressor_pop - p_expressor_ceu) × 100
```

A **positive** `delta_vs_ceu_pp` means the population has more CYP3A5 expressors than the European CEU baseline. More expressors implies more patients with higher clearance and greater divergence from European-calibrated standard dosing protocols.

The CEU baseline is `1kg:ceu` in gnomAD r3. All computations use the `gnomad_r3` dataset.

---

## Project Structure

```
latam-equity-dashboard/
├── ingest/
│   ├── _common.py               # HTTP utilities, retry logic, file I/O helpers
│   ├── gnomad_pgx.py            # Fetch CYP3A5*3 from gnomAD v3 via GraphQL
│   ├── pgx_artifacts.py         # Fetch CPIC drug impact summary (local or GitHub)
│   └── transplant_artifacts.py  # Fetch IRODaT world transplants (local or GitHub)
├── src/latam_equity/
│   ├── equity_index.py          # 4-country index (1000G + CPIC composite gap)
│   ├── equity_index_extended.py # 18-country index (gnomAD HWE, AMR proxy)
│   ├── pgx_composite.py         # Composite PGx gap across 5 CPIC Level A pairs
│   ├── timeseries.py            # Per-country transplant time series builder
│   └── export.py                # Orchestrates full pipeline → JSON artifacts
├── tests/
│   ├── test_equity_index.py
│   ├── test_equity_index_extended.py
│   └── test_timeseries.py
├── schemas/
│   └── latam-equity-index.schema.json  # JSON Schema for the primary artifact
├── dbt_project/                         # Optional dbt/DuckDB transformation layer
│   └── models/
│       ├── staging/
│       └── marts/fct_latam_equity_index.sql
├── data/
│   └── raw/        # Downloaded source files (gitignored; populated by make ingest)
├── artifacts/      # Pipeline outputs (gitignored; populated by make export)
├── Makefile
└── pyproject.toml
```

---

## Installation

Requires Python 3.11 or later.

```bash
# Clone the repository
git clone https://github.com/enriqew/latam-equity-dashboard.git
cd latam-equity-dashboard

# Install in editable mode with development dependencies
pip install -e ".[dev]"
```

---

## Running the Pipeline

The full pipeline runs in three sequential steps, each wrapped by a Makefile target.

### Step 1 — Ingest raw data

Downloads or copies source files into `data/raw/`.

```bash
make ingest
```

This is equivalent to running three sub-steps in order:

```bash
make ingest-transplants   # → data/raw/world-transplants.json
make ingest-pgx           # → data/raw/drug_impact_summary.json
                          #   data/raw/actionability_ranking.json
make ingest-gnomad        # → data/raw/gnomad_cyp3a5.json  (live gnomAD API call)
```

By default, `ingest-transplants` and `ingest-pgx` read from a local sibling repository. To fetch from GitHub instead:

```bash
python -m ingest.transplant_artifacts --source github
python -m ingest.pgx_artifacts --source github
```

### Step 2 — Export artifacts

Reads from `data/raw/`, computes both indexes and the time series, and writes JSON to `artifacts/`.

```bash
make export
```

### Step 3 — Run tests

```bash
make test
```

### Optional — dbt transformation layer

A dbt project is included for users who want to run the index computation inside DuckDB rather than in Python. It is not required to produce the JSON artifacts.

```bash
make dbt-run
make dbt-test
```

### Full pipeline (ingest + dbt + export)

```bash
make pipeline
```

### Code quality

```bash
make lint   # ruff check + mypy
```

---

## Artifact Schema

All artifacts are written to `artifacts/` as UTF-8 JSON. The directory is gitignored; files are populated at runtime by `make export`.

### `latam_equity_index.json`

Array of records, one per country (4 records). PGx gap is the composite mean across 5 CPIC Level A pairs. Sorted descending by `equity_score`.

Key fields:

| Field | Type | Description |
|-------|------|-------------|
| `population_code` | string | 1000 Genomes code: MXL, PEL, CLM, or PUR |
| `country_iso3` | string | ISO 3166-1 alpha-3 |
| `country_name` | string | Human-readable country name |
| `transplants_pmp` | number or null | Kidney transplants per million population, most recent year |
| `total_transplants` | integer or null | Absolute transplant count for the reported year |
| `data_year` | integer or null | Year of the transplant data used |
| `transplants_pmp_is_estimate` | boolean | True only for PRI (no IRODaT data available) |
| `pgx_delta_tacrolimus` | number | `delta_vs_baseline` (pp) for tacrolimus/CYP3A5 specifically |
| `pct_requiring_change` | number | % of population requiring dose/therapy change (tacrolimus/CYP3A5) |
| `pgx_composite_gap` | number | Mean `|delta_vs_baseline|` across all 5 focal pairs |
| `pgx_per_gene` | object | `{"GENE": abs_delta, ...}` breakdown for each focal pair |
| `pgx_gap` | number | Equal to `pgx_composite_gap`; named for formula symmetry |
| `access_gap` | number | `global_median_pmp - transplants_pmp`; positive = below median |
| `equity_score` | number | `access_gap + pgx_gap`; higher = worse compounded inequity |

Full JSON Schema: `schemas/latam-equity-index.schema.json`

### `latam_equity_index_extended.json`

Array of records, one per country (18 records). PGx gap is `|delta_vs_ceu_pp|` for CYP3A5\*3 under Hardy-Weinberg. Sorted descending by `equity_score`.

Fields in common with the original index: `population_code`, `country_iso3`, `country_name`, `transplants_pmp`, `total_transplants`, `data_year`, `transplants_pmp_is_estimate`, `access_gap`, `pgx_gap`, `equity_score`.

Additional fields:

| Field | Type | Description |
|-------|------|-------------|
| `pct_expressor` | number | % of population with at least one functional CYP3A5\*1 allele (HWE) |
| `delta_vs_ceu_pp` | number | `pct_expressor_pop - pct_expressor_ceu` in percentage points |
| `pgx_data_source` | string | `"gnomad_r3_hwe"` for 1000G populations; `"gnomad_r3_hwe_amr_proxy"` for AMR-proxy countries |
| `pgx_is_proxy` | boolean | True for the 14 countries using the gnomAD AMR aggregate |

### `latam_transplant_timeseries.json`

Array of country objects, sorted by `country_iso3`. Includes LATAM countries (MEX, PER, COL) and European reference countries (ESP, GBR) for visual benchmarking. Reference countries are flagged with `is_reference: true`.

```json
{
  "country_iso3": "MEX",
  "country_name": "Mexico",
  "population_code": "MXL",
  "is_reference": false,
  "series": [
    {
      "year": 2019,
      "transplants_pmp": 14.3,
      "total_transplants": 1812,
      "deceased_donors_pmp": 2.1,
      "living_donors_pmp": 9.8
    }
  ]
}
```

Rows with `transplants_pmp: null` are excluded. Series are sorted ascending by year.

### `latam_equity_meta.json`

Run metadata: generation timestamp (UTC ISO 8601), source file paths relative to the repo root, and record counts.

---

## Running Tests

The test suite uses pytest and covers 27 test cases across three modules. Tests are self-contained with synthetic fixture data and do not require network access or pre-ingested files.

```bash
make test
```

For a full coverage report:

```bash
pytest tests/ --cov=src --cov=ingest --cov-report=term-missing
```

---

## Clinical Disclaimer

The equity scores and pharmacogenomic gaps in this dataset are **epidemiological, population-level estimates**. They are computed from aggregate allele frequencies and registry-level transplant volumes. They do not represent clinical guidance, diagnostic results, or treatment recommendations for any individual patient.

CYP3A5 expressor probability derived from Hardy-Weinberg equilibrium is a statistical expectation across a population, not a genotype test result. Actual dosing decisions must be based on individual genotyping and interpreted by qualified clinicians following current CPIC guidelines.

The gnomAD AMR aggregate used as a pharmacogenomic proxy for 14 countries encompasses substantial heterogeneity across national populations. Records with `pgx_is_proxy: true` should be interpreted with that limitation in mind.

This project is not affiliated with IRODaT, the Global Observatory on Donation and Transplantation, CPIC, or the Broad Institute.

---

## References

- IRODaT / GODT: Global Observatory on Donation and Transplantation — https://www.irodat.org
- CPIC: Clinical Pharmacogenomics Implementation Consortium — https://cpicpgx.org
- 1000 Genomes Project Phase 3 — https://www.internationalgenome.org
- gnomAD v3: Genome Aggregation Database — https://gnomad.broadinstitute.org
