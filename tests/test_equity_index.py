"""Unit tests for equity index computation."""
from src.latam_equity.equity_index import build_equity_points, EquityPoint, POPULATION_TO_ISO3
from src.latam_equity.pgx_composite import CPIC_FOCAL_PAIRS


SAMPLE_TRANSPLANTS = [
    {"country_iso3": "MEX", "country_name": "Mexico",         "year": 2023, "transplants_pmp": 3.2,  "total_transplants": 400,  "source": "IRODaT"},
    {"country_iso3": "PER", "country_name": "Peru",           "year": 2023, "transplants_pmp": 0.8,  "total_transplants": 26,   "source": "IRODaT"},
    {"country_iso3": "COL", "country_name": "Colombia",       "year": 2023, "transplants_pmp": 5.1,  "total_transplants": 260,  "source": "IRODaT"},
    {"country_iso3": "ESP", "country_name": "Spain",          "year": 2023, "transplants_pmp": 60.3, "total_transplants": 2800, "source": "GODT"},
    {"country_iso3": "GBR", "country_name": "United Kingdom", "year": 2023, "transplants_pmp": 48.1, "total_transplants": 3200, "source": "GODT"},
]

SAMPLE_PGX = [
    # tacrolimus / CYP3A5
    {"drug_name": "tacrolimus",   "gene_symbol": "CYP3A5", "population_code": "MXL", "delta_vs_baseline": -2.84, "percentage_requiring_change": 6.25},
    {"drug_name": "tacrolimus",   "gene_symbol": "CYP3A5", "population_code": "PEL", "delta_vs_baseline": -6.74, "percentage_requiring_change": 2.35},
    {"drug_name": "tacrolimus",   "gene_symbol": "CYP3A5", "population_code": "CLM", "delta_vs_baseline": -4.83, "percentage_requiring_change": 4.26},
    {"drug_name": "tacrolimus",   "gene_symbol": "CYP3A5", "population_code": "CEU", "delta_vs_baseline":  0.00, "percentage_requiring_change": 9.09},
    # clopidogrel / CYP2C19
    {"drug_name": "clopidogrel",  "gene_symbol": "CYP2C19", "population_code": "MXL", "delta_vs_baseline": -21.10, "percentage_requiring_change": 20.31},
    {"drug_name": "clopidogrel",  "gene_symbol": "CYP2C19", "population_code": "PEL", "delta_vs_baseline": -33.17, "percentage_requiring_change": 8.24},
    {"drug_name": "clopidogrel",  "gene_symbol": "CYP2C19", "population_code": "CLM", "delta_vs_baseline": -16.94, "percentage_requiring_change": 24.47},
    {"drug_name": "clopidogrel",  "gene_symbol": "CYP2C19", "population_code": "CEU", "delta_vs_baseline":   0.00, "percentage_requiring_change": 41.41},
    # warfarin / VKORC1
    {"drug_name": "warfarin",     "gene_symbol": "VKORC1",  "population_code": "MXL", "delta_vs_baseline":  4.27, "percentage_requiring_change": 57.81},
    {"drug_name": "warfarin",     "gene_symbol": "VKORC1",  "population_code": "PEL", "delta_vs_baseline": 27.64, "percentage_requiring_change": 81.18},
    {"drug_name": "warfarin",     "gene_symbol": "VKORC1",  "population_code": "CLM", "delta_vs_baseline":  4.97, "percentage_requiring_change": 58.51},
    {"drug_name": "warfarin",     "gene_symbol": "VKORC1",  "population_code": "CEU", "delta_vs_baseline":  0.00, "percentage_requiring_change": 53.54},
    # fluorouracil / DPYD
    {"drug_name": "fluorouracil", "gene_symbol": "DPYD",    "population_code": "MXL", "delta_vs_baseline":  8.76, "percentage_requiring_change": 39.06},
    {"drug_name": "fluorouracil", "gene_symbol": "DPYD",    "population_code": "PEL", "delta_vs_baseline": 16.76, "percentage_requiring_change": 47.06},
    {"drug_name": "fluorouracil", "gene_symbol": "DPYD",    "population_code": "CLM", "delta_vs_baseline":  2.68, "percentage_requiring_change": 32.98},
    {"drug_name": "fluorouracil", "gene_symbol": "DPYD",    "population_code": "CEU", "delta_vs_baseline":  0.00, "percentage_requiring_change": 30.30},
    # fluvastatin / CYP2C9
    {"drug_name": "fluvastatin",  "gene_symbol": "CYP2C9",  "population_code": "MXL", "delta_vs_baseline": -8.44, "percentage_requiring_change":  4.69},
    {"drug_name": "fluvastatin",  "gene_symbol": "CYP2C9",  "population_code": "PEL", "delta_vs_baseline":-10.78, "percentage_requiring_change":  2.35},
    {"drug_name": "fluvastatin",  "gene_symbol": "CYP2C9",  "population_code": "CLM", "delta_vs_baseline": -0.36, "percentage_requiring_change": 12.77},
    {"drug_name": "fluvastatin",  "gene_symbol": "CYP2C9",  "population_code": "CEU", "delta_vs_baseline":  0.00, "percentage_requiring_change": 13.13},
]


def test_all_latam_pops_covered():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    pops = {p.population_code for p in points}
    assert pops >= {"MXL", "PEL", "CLM"}


def test_equity_score_finite():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        assert isinstance(p.equity_score, float)
        assert p.equity_score == p.equity_score  # not NaN


def test_population_map_completeness():
    assert set(POPULATION_TO_ISO3.keys()) == {"MXL", "PEL", "CLM", "PUR"}


def test_sorted_descending():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    scores = [p.equity_score for p in points]
    assert scores == sorted(scores, reverse=True)


def test_access_gap_sign():
    """Points below global median must have positive access_gap."""
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    # Global median of [3.2, 0.8, 5.1, 60.3, 48.1] = 5.1
    # MEX (3.2) is below median → access_gap > 0
    mex = next(p for p in points if p.population_code == "MXL")
    assert mex.access_gap > 0


def test_pgx_gap_is_absolute():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        assert p.pgx_gap >= 0.0


def test_pgx_gap_equals_composite():
    """pgx_gap must equal pgx_composite_gap."""
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        assert p.pgx_gap == p.pgx_composite_gap


def test_composite_gap_uses_all_focal_genes():
    """All 5 focal genes must appear in per_gene breakdown when data is present."""
    focal_genes = {gene for _, gene in CPIC_FOCAL_PAIRS}
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    mxl = next(p for p in points if p.population_code == "MXL")
    assert set(mxl.pgx_per_gene.keys()) == focal_genes


def test_composite_gap_is_mean_of_per_gene():
    """composite_gap must equal the mean of pgx_per_gene values."""
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        if p.pgx_per_gene:
            expected = round(sum(p.pgx_per_gene.values()) / len(p.pgx_per_gene), 4)
            assert p.pgx_composite_gap == expected


def test_missing_transplant_data():
    from src.latam_equity.equity_index import PRI_TRANSPLANTS_PMP_ESTIMATE
    points = build_equity_points([], SAMPLE_PGX)
    assert len(points) == len(POPULATION_TO_ISO3)
    for p in points:
        if p.country_iso3 == "PRI":
            assert p.transplants_pmp == PRI_TRANSPLANTS_PMP_ESTIMATE
            assert p.transplants_pmp_is_estimate is True
        else:
            assert p.transplants_pmp is None
            assert p.transplants_pmp_is_estimate is False
