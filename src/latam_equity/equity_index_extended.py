"""Extended LATAM Equity Index using gnomAD v3 HWE methodology for all 8 countries.

All PGx metrics are computed consistently via Hardy-Weinberg from gnomAD r3
CYP3A5*3 allele frequencies (see ingest/gnomad_pgx.py).  This allows direct
comparison across all countries on the same methodological scale.

Population → Country mapping:
  1000G pops (gnomAD 1kg:*): MXL→MEX, PEL→PER, CLM→COL, PUR→PRI
  gnomAD AMR (proxy):         AMR→ARG, BRA, CHL, URY
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from src.latam_equity.equity_index import PRI_TRANSPLANTS_PMP_ESTIMATE

# 1000G populations with direct country mapping
LATAM_1KG: dict[str, str] = {
    "1kg:mxl": "MEX",
    "1kg:pel": "PER",
    "1kg:clm": "COL",
    "1kg:pur": "PRI",
}

# New countries using gnomAD AMR as PGx proxy
GNOMAD_AMR_COUNTRIES: list[str] = ["ARG", "BRA", "CHL", "URY"]

# Display population code for AMR-proxy countries
AMR_PROXY_POP_CODE = "AMR_gnomAD"


@dataclass
class EquityPointExtended:
    population_code: str
    country_iso3: str
    country_name: str
    transplants_pmp: float | None
    total_transplants: int | None
    data_year: int | None
    pct_expressor: float          # % with ≥1 functional CYP3A5*1 (gnomAD HWE)
    delta_vs_ceu_pp: float        # delta from CEU baseline in pp
    pgx_data_source: str          # "gnomad_r3_hwe"
    pgx_is_proxy: bool            # True = AMR proxy, not direct population match
    transplants_pmp_is_estimate: bool = False
    access_gap: float = field(init=False)
    pgx_gap: float = field(init=False)
    equity_score: float = field(init=False)

    def compute(self, global_median_pmp: float) -> None:
        pmp = self.transplants_pmp or 0.0
        self.access_gap = global_median_pmp - pmp
        self.pgx_gap = abs(self.delta_vs_ceu_pp)
        self.equity_score = round(self.access_gap + self.pgx_gap, 3)


def _pgx_by_gnomad_id(gnomad_data: dict) -> dict[str, dict]:
    return {row["gnomad_pop_id"]: row for row in gnomad_data.get("populations", [])}


def build_equity_points_extended(
    world_transplants: list[dict],
    gnomad_cyp3a5: dict,
) -> list[EquityPointExtended]:
    """
    Build equity points for all 8 LATAM countries using gnomAD HWE PGx metrics.

    Args:
        world_transplants: rows from world-transplants.json
        gnomad_cyp3a5: parsed gnomad_cyp3a5.json (from ingest/gnomad_pgx.py)
    Returns:
        EquityPointExtended list sorted by equity_score DESC
    """
    # Most recent transplants_pmp per country
    latest: dict[str, dict] = {}
    for row in world_transplants:
        iso3 = row.get("country_iso3", "")
        year = int(row.get("year") or 0)
        if row.get("transplants_pmp") is None:
            continue
        if iso3 not in latest or year > int(latest[iso3].get("year") or 0):
            latest[iso3] = row

    all_pmp = [float(r["transplants_pmp"]) for r in latest.values()]
    global_median = statistics.median(all_pmp) if all_pmp else 0.0

    pgx = _pgx_by_gnomad_id(gnomad_cyp3a5)
    points: list[EquityPointExtended] = []

    # 1000G-mapped countries
    for gnomad_pop_id, iso3 in LATAM_1KG.items():
        pgx_row = pgx.get(gnomad_pop_id)
        transplant_row = latest.get(iso3)
        use_pri_estimate = iso3 == "PRI" and transplant_row is None

        pt = EquityPointExtended(
            population_code=gnomad_pop_id.replace("1kg:", "").upper(),
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else iso3,
            transplants_pmp=(
                PRI_TRANSPLANTS_PMP_ESTIMATE if use_pri_estimate
                else float(transplant_row["transplants_pmp"]) if transplant_row else None
            ),
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pct_expressor=pgx_row["pct_expressor"] if pgx_row else 0.0,
            delta_vs_ceu_pp=pgx_row["delta_vs_ceu_pp"] if pgx_row else 0.0,
            pgx_data_source="gnomad_r3_hwe",
            pgx_is_proxy=False,
            transplants_pmp_is_estimate=use_pri_estimate,
        )
        pt.compute(global_median)
        points.append(pt)

    # gnomAD AMR-proxy countries
    amr_row = pgx.get("amr")
    for iso3 in GNOMAD_AMR_COUNTRIES:
        transplant_row = latest.get(iso3)
        pt = EquityPointExtended(
            population_code=AMR_PROXY_POP_CODE,
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else iso3,
            transplants_pmp=float(transplant_row["transplants_pmp"]) if transplant_row else None,
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pct_expressor=amr_row["pct_expressor"] if amr_row else 0.0,
            delta_vs_ceu_pp=amr_row["delta_vs_ceu_pp"] if amr_row else 0.0,
            pgx_data_source="gnomad_r3_hwe_amr_proxy",
            pgx_is_proxy=True,
            transplants_pmp_is_estimate=False,
        )
        pt.compute(global_median)
        points.append(pt)

    return sorted(points, key=lambda p: p.equity_score, reverse=True)
