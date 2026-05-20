"""Unit tests for equity index computation."""
import pytest
from src.latam_equity.equity_index import build_equity_points, EquityPoint, POPULATION_TO_ISO3


SAMPLE_TRANSPLANTS = [
    {"country_iso3": "MEX", "country_name": "Mexico", "year": 2023, "transplants_pmp": 3.2, "total_transplants": 400, "source": "IRODaT"},
    {"country_iso3": "PER", "country_name": "Peru", "year": 2023, "transplants_pmp": 0.8, "total_transplants": 26, "source": "IRODaT"},
    {"country_iso3": "COL", "country_name": "Colombia", "year": 2023, "transplants_pmp": 5.1, "total_transplants": 260, "source": "IRODaT"},
    {"country_iso3": "ESP", "country_name": "Spain", "year": 2023, "transplants_pmp": 60.3, "total_transplants": 2800, "source": "GODT"},
    {"country_iso3": "GBR", "country_name": "United Kingdom", "year": 2023, "transplants_pmp": 48.1, "total_transplants": 3200, "source": "GODT"},
]

SAMPLE_PGX = [
    {"drug_name": "tacrolimus", "gene_symbol": "CYP3A5", "population_code": "MXL", "delta_vs_baseline": -2.84, "percentage_requiring_change": 6.25},
    {"drug_name": "tacrolimus", "gene_symbol": "CYP3A5", "population_code": "PEL", "delta_vs_baseline": -6.74, "percentage_requiring_change": 2.35},
    {"drug_name": "tacrolimus", "gene_symbol": "CYP3A5", "population_code": "CLM", "delta_vs_baseline": -4.83, "percentage_requiring_change": 4.26},
    {"drug_name": "tacrolimus", "gene_symbol": "CYP3A5", "population_code": "CEU", "delta_vs_baseline": 0.0, "percentage_requiring_change": 9.09},
]


def test_all_latam_pops_covered():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    pops = {p.population_code for p in points}
    assert pops >= {"MXL", "PEL", "CLM"}  # PUR/PRI may be absent from sample


def test_equity_score_finite():
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        assert isinstance(p.equity_score, float)
        assert not (p.equity_score != p.equity_score)  # not NaN


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
    """pgx_gap must be non-negative regardless of delta sign."""
    points = build_equity_points(SAMPLE_TRANSPLANTS, SAMPLE_PGX)
    for p in points:
        assert p.pgx_gap >= 0.0


def test_missing_transplant_data():
    """A population with no transplant match should still produce an EquityPoint."""
    points = build_equity_points([], SAMPLE_PGX)
    # No transplant data → all populations produce points with transplants_pmp=None
    assert len(points) == len(POPULATION_TO_ISO3)
    for p in points:
        assert p.transplants_pmp is None
