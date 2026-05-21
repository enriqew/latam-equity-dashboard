"""Unit tests for transplant time series builder."""
from src.latam_equity.timeseries import build_timeseries, LATAM_COUNTRIES, REFERENCE_COUNTRIES

SAMPLE_TRANSPLANTS = [
    {"country_iso3": "MEX", "country_name": "Mexico", "year": 2020, "transplants_pmp": 7.9,  "total_transplants": 998,  "deceased_donors_pmp": 1.4, "living_donors_pmp": 0.0},
    {"country_iso3": "MEX", "country_name": "Mexico", "year": 2021, "transplants_pmp": 16.8, "total_transplants": 2141, "deceased_donors_pmp": 2.0, "living_donors_pmp": 11.8},
    {"country_iso3": "PER", "country_name": "Peru",   "year": 2021, "transplants_pmp": 4.1,  "total_transplants": 134,  "deceased_donors_pmp": 0.8, "living_donors_pmp": 3.3},
    {"country_iso3": "COL", "country_name": "Colombia","year": 2021, "transplants_pmp": 20.0, "total_transplants": 1000, "deceased_donors_pmp": 8.0, "living_donors_pmp": 12.0},
    {"country_iso3": "ESP", "country_name": "Spain",  "year": 2021, "transplants_pmp": 60.0, "total_transplants": 2800, "deceased_donors_pmp": 40.0,"living_donors_pmp": 20.0},
    {"country_iso3": "GBR", "country_name": "United Kingdom", "year": 2021, "transplants_pmp": 48.0, "total_transplants": 3200, "deceased_donors_pmp": 30.0, "living_donors_pmp": 18.0},
    # null pmp row — should be excluded
    {"country_iso3": "MEX", "country_name": "Mexico", "year": 2019, "transplants_pmp": None, "total_transplants": 800, "deceased_donors_pmp": None, "living_donors_pmp": None},
    # unrelated country — should be excluded
    {"country_iso3": "USA", "country_name": "United States", "year": 2021, "transplants_pmp": 90.0, "total_transplants": 30000, "deceased_donors_pmp": 60.0, "living_donors_pmp": 30.0},
]


def test_returns_all_target_countries():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    iso3s = {e["country_iso3"] for e in result}
    assert iso3s == set(LATAM_COUNTRIES) | REFERENCE_COUNTRIES


def test_excludes_null_pmp_rows():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    mex = next(e for e in result if e["country_iso3"] == "MEX")
    years = [p["year"] for p in mex["series"]]
    assert 2019 not in years  # the null-pmp row must be excluded


def test_excludes_unrelated_countries():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    assert all(e["country_iso3"] != "USA" for e in result)


def test_series_sorted_by_year():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    for entry in result:
        years = [p["year"] for p in entry["series"]]
        assert years == sorted(years)


def test_latam_has_population_code():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    for entry in result:
        if entry["country_iso3"] in LATAM_COUNTRIES:
            assert entry["population_code"] == LATAM_COUNTRIES[entry["country_iso3"]]
            assert entry["is_reference"] is False


def test_reference_countries_flagged():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    for entry in result:
        if entry["country_iso3"] in REFERENCE_COUNTRIES:
            assert entry["population_code"] is None
            assert entry["is_reference"] is True


def test_series_contains_expected_fields():
    result = build_timeseries(SAMPLE_TRANSPLANTS)
    mex = next(e for e in result if e["country_iso3"] == "MEX")
    for point in mex["series"]:
        assert "year" in point
        assert "transplants_pmp" in point
        assert "deceased_donors_pmp" in point
        assert "living_donors_pmp" in point


def test_empty_input_returns_empty():
    assert build_timeseries([]) == []
