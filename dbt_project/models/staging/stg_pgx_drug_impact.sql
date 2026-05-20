{{ config(materialized='view') }}

/*
  stg_pgx_drug_impact
  ────────────────────
  Source: data/raw/drug_impact_summary.json  (pgx-latam-atlas pipeline output)

  Columns:
    drug_name                 — e.g. "tacrolimus", "warfarin"
    gene_symbol               — e.g. "CYP3A5", "VKORC1"
    population_code           — 1000 Genomes code: MXL, PEL, CLM, PUR, CEU …
    percentage_requiring_change — % of individuals needing dose/therapy change
    baseline_ceu_percentage   — equivalent percentage in CEU (European proxy)
    delta_vs_baseline         — percentage_requiring_change − baseline_ceu_percentage (pp)
*/

WITH raw AS (
    SELECT *
    FROM read_json_auto(
        '{{ env_var("REPO_ROOT", "..") }}/data/raw/drug_impact_summary.json',
        auto_detect = true
    )
)

SELECT
    CAST(drug_name                    AS VARCHAR) AS drug_name,
    CAST(gene_symbol                  AS VARCHAR) AS gene_symbol,
    CAST(population_code              AS VARCHAR) AS population_code,
    TRY_CAST(percentage_requiring_change AS DOUBLE) AS percentage_requiring_change,
    TRY_CAST(baseline_ceu_percentage  AS DOUBLE)  AS baseline_ceu_percentage,
    TRY_CAST(delta_vs_baseline        AS DOUBLE)  AS delta_vs_baseline
FROM raw
WHERE drug_name       IS NOT NULL
  AND gene_symbol     IS NOT NULL
  AND population_code IS NOT NULL
