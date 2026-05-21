"""Tests for the extended equity index (gnomAD multi-gene HWE methodology)."""
import dataclasses
import pytest
from src.latam_equity.equity_index_extended import (
    build_equity_points_extended,
    LATAM_1KG,
    GNOMAD_AMR_COUNTRIES,
    AMR_PROXY_POP_CODE,
)

SAMPLE_TRANSPLANTS = [
    # 1000G-mapped countries
    {"country_iso3": "MEX", "country_name": "Mexico",      "year": 2024, "transplants_pmp": 23.5, "total_transplants": 3080},
    {"country_iso3": "PER", "country_name": "Peru",        "year": 2024, "transplants_pmp": 6.3,  "total_transplants": 217},
    {"country_iso3": "COL", "country_name": "Colombia",    "year": 2024, "transplants_pmp": 26.3, "total_transplants": 1391},
    # AMR-proxy countries (original set)
    {"country_iso3": "ARG", "country_name": "Argentina",   "year": 2024, "transplants_pmp": 48.8, "total_transplants": 2200},
    {"country_iso3": "BRA", "country_name": "Brazil",      "year": 2024, "transplants_pmp": 42.0, "total_transplants": 8900},
    {"country_iso3": "CHL", "country_name": "Chile",       "year": 2024, "transplants_pmp": 35.6, "total_transplants": 650},
    {"country_iso3": "URY", "country_name": "Uruguay",     "year": 2024, "transplants_pmp": 49.6, "total_transplants": 175},
    # AMR-proxy countries (expanded)
    {"country_iso3": "BOL", "country_name": "Bolivia",     "year": 2024, "transplants_pmp": 7.8,  "total_transplants": 95},
    {"country_iso3": "CRI", "country_name": "Costa Rica",  "year": 2024, "transplants_pmp": 57.0, "total_transplants": 290},
    {"country_iso3": "ECU", "country_name": "Ecuador",     "year": 2024, "transplants_pmp": 19.3, "total_transplants": 350},
    {"country_iso3": "GTM", "country_name": "Guatemala",   "year": 2024, "transplants_pmp": 15.0, "total_transplants": 240},
    {"country_iso3": "HND", "country_name": "Honduras",    "year": 2024, "transplants_pmp": 0.8,  "total_transplants": 10},
    {"country_iso3": "NIC", "country_name": "Nicaragua",   "year": 2024, "transplants_pmp": 6.3,  "total_transplants": 42},
    {"country_iso3": "PAN", "country_name": "Panama",      "year": 2024, "transplants_pmp": 19.2, "total_transplants": 80},
    {"country_iso3": "PRY", "country_name": "Paraguay",    "year": 2024, "transplants_pmp": 8.8,  "total_transplants": 65},
    {"country_iso3": "SLV", "country_name": "El Salvador", "year": 2024, "transplants_pmp": 8.0,  "total_transplants": 50},
    {"country_iso3": "VEN", "country_name": "Venezuela",   "year": 2023, "transplants_pmp": 15.1, "total_transplants": 450},
    # Reference (not in index, but affects global_median)
    {"country_iso3": "ESP", "country_name": "Spain",       "year": 2024, "transplants_pmp": 60.3, "total_transplants": 2800},
]

# gnomad_pgx_multi.json format — 5 genes, 6 populations each
SAMPLE_GNOMAD_MULTI = {
    "generated_at": "2025-01-01T00:00:00+00:00",
    "dataset": "gnomad_r3",
    "methodology": "HWE p_lof_carrier = 1-(1-AF_lof)^2.",
    "genes": {
        "CYP3A5": {
            "variant_id": "7-99672916-T-C", "rsid": "rs776746",
            "drug": "tacrolimus", "star_allele": "*3",
            "populations": [
                {"gnomad_pop_id": "1kg:ceu", "af_lof": 0.9622, "p_lof_carrier": 0.9986, "pct_lof_carrier": 99.86, "delta_vs_ceu_pp": 0.0},
                {"gnomad_pop_id": "1kg:mxl", "af_lof": 0.7540, "p_lof_carrier": 0.9395, "pct_lof_carrier": 93.95, "delta_vs_ceu_pp": -5.91},
                {"gnomad_pop_id": "1kg:clm", "af_lof": 0.8263, "p_lof_carrier": 0.9698, "pct_lof_carrier": 96.98, "delta_vs_ceu_pp": -2.88},
                {"gnomad_pop_id": "1kg:pel", "af_lof": 0.8779, "p_lof_carrier": 0.9851, "pct_lof_carrier": 98.51, "delta_vs_ceu_pp": -1.35},
                {"gnomad_pop_id": "1kg:pur", "af_lof": 0.7374, "p_lof_carrier": 0.9310, "pct_lof_carrier": 93.10, "delta_vs_ceu_pp": -6.76},
                {"gnomad_pop_id": "amr",     "af_lof": 0.7795, "p_lof_carrier": 0.9514, "pct_lof_carrier": 95.14, "delta_vs_ceu_pp": -4.72},
            ],
        },
        "CYP2C19": {
            "variant_id": "10-94781859-G-A", "rsid": "rs4244285",
            "drug": "clopidogrel", "star_allele": "*2",
            "populations": [
                {"gnomad_pop_id": "1kg:ceu", "af_lof": 0.1303, "p_lof_carrier": 0.2435, "pct_lof_carrier": 24.35, "delta_vs_ceu_pp": 0.0},
                {"gnomad_pop_id": "1kg:mxl", "af_lof": 0.1270, "p_lof_carrier": 0.2378, "pct_lof_carrier": 23.78, "delta_vs_ceu_pp": -0.57},
                {"gnomad_pop_id": "1kg:clm", "af_lof": 0.1053, "p_lof_carrier": 0.1994, "pct_lof_carrier": 19.94, "delta_vs_ceu_pp": -4.41},
                {"gnomad_pop_id": "1kg:pel", "af_lof": 0.0640, "p_lof_carrier": 0.1238, "pct_lof_carrier": 12.38, "delta_vs_ceu_pp": -11.97},
                {"gnomad_pop_id": "1kg:pur", "af_lof": 0.1313, "p_lof_carrier": 0.2454, "pct_lof_carrier": 24.54, "delta_vs_ceu_pp": 0.19},
                {"gnomad_pop_id": "amr",     "af_lof": 0.1302, "p_lof_carrier": 0.2435, "pct_lof_carrier": 24.35, "delta_vs_ceu_pp": 0.00},
            ],
        },
        "VKORC1": {
            "variant_id": "16-31096368-C-T", "rsid": "rs9923231",
            "drug": "warfarin", "star_allele": "-1639A",
            "populations": [
                {"gnomad_pop_id": "1kg:ceu", "af_lof": 0.4034, "p_lof_carrier": 0.6440, "pct_lof_carrier": 64.40, "delta_vs_ceu_pp": 0.0},
                {"gnomad_pop_id": "1kg:mxl", "af_lof": 0.4683, "p_lof_carrier": 0.7172, "pct_lof_carrier": 71.72, "delta_vs_ceu_pp": 7.32},
                {"gnomad_pop_id": "1kg:clm", "af_lof": 0.4158, "p_lof_carrier": 0.6587, "pct_lof_carrier": 65.87, "delta_vs_ceu_pp": 1.47},
                {"gnomad_pop_id": "1kg:pel", "af_lof": 0.3837, "p_lof_carrier": 0.6202, "pct_lof_carrier": 62.02, "delta_vs_ceu_pp": -2.38},
                {"gnomad_pop_id": "1kg:pur", "af_lof": 0.3788, "p_lof_carrier": 0.6141, "pct_lof_carrier": 61.41, "delta_vs_ceu_pp": -2.99},
                {"gnomad_pop_id": "amr",     "af_lof": 0.3910, "p_lof_carrier": 0.6291, "pct_lof_carrier": 62.91, "delta_vs_ceu_pp": -1.49},
            ],
        },
        "DPYD": {
            "variant_id": "1-97450058-C-T", "rsid": "rs3918290",
            "drug": "fluorouracil", "star_allele": "*2A",
            "populations": [
                {"gnomad_pop_id": "1kg:ceu", "af_lof": 0.0042, "p_lof_carrier": 0.0084, "pct_lof_carrier": 0.84, "delta_vs_ceu_pp": 0.0},
                {"gnomad_pop_id": "1kg:mxl", "af_lof": 0.0000, "p_lof_carrier": 0.0000, "pct_lof_carrier": 0.00, "delta_vs_ceu_pp": -0.84},
                {"gnomad_pop_id": "1kg:clm", "af_lof": 0.0000, "p_lof_carrier": 0.0000, "pct_lof_carrier": 0.00, "delta_vs_ceu_pp": -0.84},
                {"gnomad_pop_id": "1kg:pel", "af_lof": 0.0058, "p_lof_carrier": 0.0116, "pct_lof_carrier": 1.16, "delta_vs_ceu_pp": 0.32},
                {"gnomad_pop_id": "1kg:pur", "af_lof": 0.0000, "p_lof_carrier": 0.0000, "pct_lof_carrier": 0.00, "delta_vs_ceu_pp": -0.84},
                {"gnomad_pop_id": "amr",     "af_lof": 0.0017, "p_lof_carrier": 0.0034, "pct_lof_carrier": 0.34, "delta_vs_ceu_pp": -0.50},
            ],
        },
        "CYP2C9": {
            "variant_id": "10-94981296-A-C", "rsid": "rs1057910",
            "drug": "fluvastatin", "star_allele": "*3",
            "populations": [
                {"gnomad_pop_id": "1kg:ceu", "af_lof": 0.0630, "p_lof_carrier": 0.1221, "pct_lof_carrier": 12.21, "delta_vs_ceu_pp": 0.0},
                {"gnomad_pop_id": "1kg:mxl", "af_lof": 0.0317, "p_lof_carrier": 0.0625, "pct_lof_carrier": 6.25,  "delta_vs_ceu_pp": -5.96},
                {"gnomad_pop_id": "1kg:clm", "af_lof": 0.0632, "p_lof_carrier": 0.1223, "pct_lof_carrier": 12.23, "delta_vs_ceu_pp": 0.02},
                {"gnomad_pop_id": "1kg:pel", "af_lof": 0.0116, "p_lof_carrier": 0.0231, "pct_lof_carrier": 2.31,  "delta_vs_ceu_pp": -9.90},
                {"gnomad_pop_id": "1kg:pur", "af_lof": 0.0404, "p_lof_carrier": 0.0792, "pct_lof_carrier": 7.92,  "delta_vs_ceu_pp": -4.29},
                {"gnomad_pop_id": "amr",     "af_lof": 0.0491, "p_lof_carrier": 0.0959, "pct_lof_carrier": 9.59,  "delta_vs_ceu_pp": -2.62},
            ],
        },
    },
}


def test_returns_all_expected_countries():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    iso3s = {p.country_iso3 for p in points}
    expected = set(LATAM_1KG.values()) | set(GNOMAD_AMR_COUNTRIES)
    assert iso3s == expected


def test_sorted_descending():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    scores = [p.equity_score for p in points]
    assert scores == sorted(scores, reverse=True)


def test_pgx_gap_non_negative():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        assert p.pgx_gap >= 0.0


def test_amr_proxy_countries_flagged():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        if p.country_iso3 in GNOMAD_AMR_COUNTRIES:
            assert p.pgx_is_proxy is True
            assert p.population_code == AMR_PROXY_POP_CODE
            assert "amr_proxy" in p.pgx_data_source


def test_1kg_countries_not_proxy():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        if p.country_iso3 in LATAM_1KG.values():
            assert p.pgx_is_proxy is False


def test_amr_proxy_countries_share_pgx_values():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    amr_points = [p for p in points if p.country_iso3 in GNOMAD_AMR_COUNTRIES]
    composites = {p.pgx_composite_gap for p in amr_points}
    assert len(composites) == 1  # all share same PGx data


def test_pri_uses_estimate_when_missing():
    from src.latam_equity.equity_index import PRI_TRANSPLANTS_PMP_ESTIMATE
    points = build_equity_points_extended([], SAMPLE_GNOMAD_MULTI)
    pri = next((p for p in points if p.country_iso3 == "PRI"), None)
    assert pri is not None
    assert pri.transplants_pmp == PRI_TRANSPLANTS_PMP_ESTIMATE
    assert pri.transplants_pmp_is_estimate is True


def test_equity_score_finite():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        assert isinstance(p.equity_score, float)
        assert p.equity_score == p.equity_score  # not NaN


def test_dataclass_serializable():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        d = dataclasses.asdict(p)
        assert "equity_score" in d
        assert "pgx_is_proxy" in d
        assert "pgx_data_source" in d
        assert "pgx_composite_gap" in d
        assert "pgx_per_gene" in d


def test_per_gene_has_all_5_genes():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    mxl = next(p for p in points if p.country_iso3 == "MEX")
    assert set(mxl.pgx_per_gene.keys()) == {"CYP3A5", "CYP2C19", "VKORC1", "DPYD", "CYP2C9"}


def test_composite_equals_mean_of_per_gene():
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    for p in points:
        if p.pgx_per_gene:
            expected = round(sum(p.pgx_per_gene.values()) / len(p.pgx_per_gene), 4)
            assert p.pgx_composite_gap == expected


def test_cyp3a5_uses_expressor_framing():
    """CYP3A5 gap should reflect expressor delta (~35 pp for MXL), not LOF-carrier delta (~6 pp)."""
    points = build_equity_points_extended(SAMPLE_TRANSPLANTS, SAMPLE_GNOMAD_MULTI)
    mxl = next(p for p in points if p.country_iso3 == "MEX")
    # Expressor delta ≈ 35.73 pp; LOF-carrier delta ≈ 5.91 pp — must use expressor framing
    assert mxl.pgx_per_gene["CYP3A5"] > 30.0
