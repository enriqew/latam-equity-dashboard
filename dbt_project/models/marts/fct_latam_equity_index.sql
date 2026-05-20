{{ config(materialized='table') }}

WITH latest_transplants AS (
    SELECT
        country_iso3,
        country_name,
        transplants_pmp,
        total_transplants,
        year,
        ROW_NUMBER() OVER (PARTITION BY country_iso3 ORDER BY year DESC) AS rn
    FROM {{ ref('stg_world_transplants') }}
    WHERE transplants_pmp IS NOT NULL
),
global_stats AS (
    SELECT
        MEDIAN(transplants_pmp) AS global_median_pmp
    FROM latest_transplants
    WHERE rn = 1
),
pgx_tacrolimus AS (
    SELECT
        population_code,
        delta_vs_baseline,
        percentage_requiring_change
    FROM {{ ref('stg_pgx_drug_impact') }}
    WHERE drug_name = 'tacrolimus'
      AND gene_symbol = 'CYP3A5'
      AND population_code != 'CEU'
),
pop_to_iso AS (
    SELECT * FROM (VALUES
        ('MXL', 'MEX'),
        ('PEL', 'PER'),
        ('CLM', 'COL'),
        ('PUR', 'PRI')
    ) AS t(population_code, country_iso3)
)
SELECT
    p.population_code,
    p.country_iso3,
    lt.country_name,
    lt.transplants_pmp,
    lt.total_transplants,
    lt.year AS data_year,
    pgx.delta_vs_baseline AS pgx_delta_tacrolimus,
    pgx.percentage_requiring_change AS pct_requiring_change,
    -- Access gap: how far below global median (positive = worse access)
    ROUND(gs.global_median_pmp - lt.transplants_pmp, 2) AS access_gap,
    -- PGx gap: |delta_vs_baseline|
    ROUND(ABS(pgx.delta_vs_baseline), 2) AS pgx_gap,
    -- Equity score: sum of both gaps (higher = worse compounded inequity)
    ROUND((gs.global_median_pmp - lt.transplants_pmp) + ABS(pgx.delta_vs_baseline), 3) AS equity_score
FROM pop_to_iso p
JOIN pgx_tacrolimus pgx ON p.population_code = pgx.population_code
LEFT JOIN latest_transplants lt ON p.country_iso3 = lt.country_iso3 AND lt.rn = 1
CROSS JOIN global_stats gs
ORDER BY equity_score DESC
