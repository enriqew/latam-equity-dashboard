{{ config(materialized='view') }}

/*
  stg_world_transplants
  ─────────────────────
  Source: data/raw/world-transplants.json  (IRODaT / GODT, 104 countries, 1993–2025)

  This staging model reads the raw JSON file via DuckDB's read_json_auto and
  casts all columns to consistent types.  Rows with NULL transplants_pmp are
  kept here; downstream models filter as needed.
*/

WITH raw AS (
    SELECT *
    FROM read_json_auto(
        '{{ env_var("REPO_ROOT", "..") }}/data/raw/world-transplants.json',
        auto_detect  = true,
        maximum_object_size = 33554432
    )
)

SELECT
    CAST(country_iso3           AS VARCHAR)  AS country_iso3,
    CAST(country_name           AS VARCHAR)  AS country_name,
    CAST(year                   AS INTEGER)  AS year,
    TRY_CAST(total_transplants  AS INTEGER)  AS total_transplants,
    TRY_CAST(transplants_pmp    AS DOUBLE)   AS transplants_pmp,
    TRY_CAST(deceased_donors_pmp AS DOUBLE)  AS deceased_donors_pmp,
    TRY_CAST(living_donors_pmp  AS DOUBLE)   AS living_donors_pmp,
    CAST(source                 AS VARCHAR)  AS source
FROM raw
WHERE country_iso3 IS NOT NULL
  AND year         IS NOT NULL
