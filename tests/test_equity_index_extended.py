"""Tests for the extended equity index (gnomAD HWE methodology)."""
import dataclasses
import pytest
from src.latam_equity.equity_index_extended import (
    build_equity_points_extended,
    LATAM_1KG,
    GNOMAD_AMR_COUNTRIES,
    AMR_PROXY_POP_CODE,
)

SAMPLE_TRANSPLANTS = [
    {"country_iso3": "MEX", "country_name": "Mexico",    "year": 2024, "transplants_pmp": 23.5, "total_transplants": 3080},
    {"country_iso3": "PER", "country_name": "Peru",      "year": 2024, "transplants_pmp": 6.3,  "total_transplants": 217},
    {"country_iso3": "COL", "country_name": "Colombia",  "year": 2024, "transplants_pmp": 26.3, "total_transplants": 1391},
    {"country_iso3": "ARG", "country_name": "Argentina", "year": 2024, "transplants_pmp": 48.8, "total_transplants": 2200},
    {"country_iso3": "BRA", "country_name": "Brazil",    "year": 2024, "transplants_pmp": 42.0, "total_transplants": 8900},
    {"country_iso3": "CHL", "country_name": "Chile",     "year": 2024, "transplants_pmp": 35.6, "total_transplants": 650},
    {"country_iso3": "URY", "country_name": "Uruguay",   "year": 2024, "transplants_pmp": 49.6, "total_transplants": 175},
    {"country_iso3": "ESP", "country_name": "Spain",     "year": 2024, "transplants_pmp": 60.3, "total_transplants": 2800},
]

SAMPLE_GNOMAD = {
    "variant_id": "7-99672916-T-C",
    "gene": "CYP3A5",
    "drug": "tacrolimus",
    "populations": [
        {"gnomad_pop_id": "1kg:ceu", "ac_star3": 229, "an": 238, "af_star3": 0.9622, "p_expressor": 0.0742, "pct_expressor": 7.42,  "delta_vs_ceu_pp": 0.0},
        {"gnomad_pop_id": "1kg:mxl", "ac_star3": 95,  "an": 126, "af_star3": 0.7540, "p_expressor": 0.4315, "pct_expressor": 43.15, "delta_vs_ceu_pp": 35.73},
        {"gnomad_pop_id": "1kg:clm", "ac_star3": 157, "an": 190, "af_star3": 0.8263, "p_expressor": 0.3172, "pct_expressor": 31.72, "delta_vs_ceu_pp": 24.30},
        {"gnomad_pop_id": "1kg:pel", "ac_star3": 151, "an": 172, "af_star3": 0.8779, "p_expressor": 0.2293, "pct_expressor": 22.93, "delta_vs_ceu_pp": 15.51},
        {"gnomad_pop_id": "1kg:pur", "ac_star3": 146, "an": 198, "af_star3": 0.7374, "p_expressor": 0.4563, "pct_expressor": 45.63, "delta_vs_ceu_pp": 38.21},
        {"gnomad_pop_id": "amr",     "ac_star3": 11892, "an": 15256, "af_star3": 0.7795, "p_expressor": 0.3924, "pct_expressor": 39.24, "delta_vs_ceu_pp": 31.82},
    ],
}


def test_returns_all_8_countries():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    iso3s = {p.country_iso3 for p in points}
    expected = set(LATAM_1KG.values()) | set(GNOMAD_AMR_COUNTRIES)
    assert iso3s == expected


def test_sorted_descending():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    scores = [p.equity_score for p in points]
    assert scores == sorted(scores, reverse=True)


def test_pgx_gap_non_negative():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    for p in points:
        assert p.pgx_gap >= 0.0


def test_amr_proxy_countries_flagged():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    for p in points:
        if p.country_iso3 in GNOMAD_AMR_COUNTRIES:
            assert p.pgx_is_proxy is True
            assert p.population_code == AMR_PROXY_POP_CODE
            assert "amr_proxy" in p.pgx_data_source


def test_1kg_countries_not_proxy():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    for p in points:
        if p.country_iso3 in LATAM_1KG.values():
            assert p.pgx_is_proxy is False


def test_amr_proxy_countries_share_pgx_values():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    amr_points = [p for p in points if p.country_iso3 in GNOMAD_AMR_COUNTRIES]
    deltas = {p.delta_vs_ceu_pp for p in amr_points}
    assert len(deltas) == 1  # all share same PGx data


def test_pri_uses_estimate_when_missing():
    from src.latam_equity.equity_index import PRI_TRANSPLANTS_PMP_ESTIMATE
    points = build_equity_points_extended([], SAMPLE_GNOMAD)
    pri = next((p for p in points if p.country_iso3 == "PRI"), None)
    assert pri is not None
    assert pri.transplants_pmp == PRI_TRANSPLANTS_PMP_ESTIMATE
    assert pri.transplants_pmp_is_estimate is True


def test_equity_score_finite():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    for p in points:
        assert isinstance(p.equity_score, float)
        assert p.equity_score == p.equity_score  # not NaN


def test_dataclass_serializable():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD)
    for p in points:
        d = dataclasses.asdict(p)
        assert "equity_score" in d
        assert "pgx_is_proxy" in d
        assert "pgx_data_source" in d
