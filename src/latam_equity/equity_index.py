"""Compute LATAM Equity Index: access gap × pharmacogenomic gap.

Equity Index components:
  - Access component: transplants_pmp relative to global median
    (negative = below median; positive = above)
  - PGx component (composite): mean abs(delta_vs_baseline) across 5 CPIC Level A
    gene-drug pairs relevant to transplant care (see pgx_composite.py)
  - Combined equity score: higher = worse compounded inequity

Population → Country mapping (1000 Genomes → ISO3):
  MXL → MEX, PEL → PER, CLM → COL, PUR → PRI
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from src.latam_equity.pgx_composite import compute_composite_gap


POPULATION_TO_ISO3: dict[str, str] = {
    "MXL": "MEX",
    "PEL": "PER",
    "CLM": "COL",
    "PUR": "PRI",
}

ISO3_TO_COUNTRY_NAME: dict[str, str] = {
    "MEX": "Mexico",
    "PER": "Peru",
    "COL": "Colombia",
    "PRI": "Puerto Rico",
}

# IRODaT does not report Puerto Rico as a separate jurisdiction.
# This estimate is derived from US territory data (~45 pmp).
PRI_TRANSPLANTS_PMP_ESTIMATE = 45.0


@dataclass
class EquityPoint:
    population_code: str
    country_iso3: str
    country_name: str
    transplants_pmp: float | None
    total_transplants: int | None
    data_year: int | None
    pgx_delta_tacrolimus: float        # delta_vs_baseline for tacrolimus/CYP3A5
    pct_requiring_change: float        # % needing dose adjustment (tacrolimus)
    pgx_composite_gap: float           # mean abs(delta) across 5 CPIC Level A pairs
    pgx_per_gene: dict[str, float]     # {gene: abs_delta} breakdown
    transplants_pmp_is_estimate: bool = False
    access_gap: float = field(init=False)
    pgx_gap: float = field(init=False)           # = pgx_composite_gap
    equity_score: float = field(init=False)      # access_gap + pgx_gap (higher = worse)

    def compute(self, global_median_pmp: float) -> None:
        pmp = self.transplants_pmp or 0.0
        self.access_gap = global_median_pmp - pmp
        self.pgx_gap = self.pgx_composite_gap
        self.equity_score = round(self.access_gap + self.pgx_gap, 3)


def build_equity_points(
    world_transplants: list[dict],
    drug_impact: list[dict],
) -> list[EquityPoint]:
    """
    Join transplant volumes with PGx delta for LATAM populations.

    Args:
        world_transplants: rows from world-transplants.json
        drug_impact: rows from drug_impact_summary.json
    Returns:
        EquityPoint list sorted by equity_score DESC
    """
    # Get most recent transplants_pmp per country
    latest: dict[str, dict] = {}
    for row in world_transplants:
        iso3 = row.get("country_iso3", "")
        year = int(row.get("year") or 0)
        if row.get("transplants_pmp") is None:
            continue
        if iso3 not in latest or year > int(latest[iso3].get("year") or 0):
            latest[iso3] = row

    # Compute global median pmp
    all_pmp = [float(r["transplants_pmp"]) for r in latest.values() if r.get("transplants_pmp")]
    global_median = statistics.median(all_pmp) if all_pmp else 0.0

    # Index tacrolimus/CYP3A5 rows by population (for backward-compat fields)
    tacro_by_pop: dict[str, dict] = {}
    for row in drug_impact:
        if row.get("drug_name") == "tacrolimus" and row.get("gene_symbol") == "CYP3A5":
            pop = row.get("population_code", "")
            if pop and pop != "CEU":
                tacro_by_pop[pop] = row

    points: list[EquityPoint] = []
    for pop, iso3 in POPULATION_TO_ISO3.items():
        tacro = tacro_by_pop.get(pop)
        transplant_row = latest.get(iso3)
        composite_gap, per_gene = compute_composite_gap(drug_impact, pop)

        use_pri_estimate = (iso3 == "PRI" and transplant_row is None)
        point = EquityPoint(
            population_code=pop,
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else ISO3_TO_COUNTRY_NAME.get(iso3, iso3),
            transplants_pmp=(
                PRI_TRANSPLANTS_PMP_ESTIMATE if use_pri_estimate
                else float(transplant_row["transplants_pmp"]) if transplant_row else None
            ),
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pgx_delta_tacrolimus=float(tacro.get("delta_vs_baseline", 0)) if tacro else 0.0,
            pct_requiring_change=float(tacro.get("percentage_requiring_change", 0)) if tacro else 0.0,
            pgx_composite_gap=composite_gap,
            pgx_per_gene=per_gene,
            transplants_pmp_is_estimate=use_pri_estimate,
        )
        point.compute(global_median)
        points.append(point)

    return sorted(points, key=lambda p: p.equity_score, reverse=True)
