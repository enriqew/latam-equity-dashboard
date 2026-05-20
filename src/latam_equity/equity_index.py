"""Compute LATAM Equity Index: access gap × pharmacogenomic gap.

Equity Index components:
  - Access component: transplants_pmp relative to global median
    (negative = below median; positive = above)
  - PGx component: |delta_vs_baseline| for tacrolimus (CYP3A5)
    (higher = larger gap from European dosing protocol)
  - Combined equity score: higher = worse compounded inequity

Population → Country mapping (1000 Genomes → ISO3):
  MXL → MEX, PEL → PER, CLM → COL, PUR → PRI
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field


POPULATION_TO_ISO3: dict[str, str] = {
    "MXL": "MEX",
    "PEL": "PER",
    "CLM": "COL",
    "PUR": "PRI",
}


@dataclass
class EquityPoint:
    population_code: str
    country_iso3: str
    country_name: str
    transplants_pmp: float | None
    total_transplants: int | None
    data_year: int | None
    pgx_delta_tacrolimus: float  # |delta_vs_baseline| CYP3A5
    pct_requiring_change: float  # % needing dose adjustment
    access_gap: float = field(init=False)      # below global median = negative
    pgx_gap: float = field(init=False)         # already == |pgx_delta_tacrolimus|
    equity_score: float = field(init=False)    # access_gap + pgx_gap (higher = worse)

    def compute(self, global_median_pmp: float) -> None:
        pmp = self.transplants_pmp or 0.0
        self.access_gap = global_median_pmp - pmp   # positive when below median
        self.pgx_gap = abs(self.pgx_delta_tacrolimus)
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

    # Get PGx delta for tacrolimus/CYP3A5 per LATAM population
    pgx_by_pop: dict[str, dict] = {}
    for row in drug_impact:
        if row.get("drug_name") == "tacrolimus" and row.get("gene_symbol") == "CYP3A5":
            pop = row.get("population_code", "")
            if pop and pop != "CEU":
                pgx_by_pop[pop] = row

    points: list[EquityPoint] = []
    for pop, iso3 in POPULATION_TO_ISO3.items():
        pgx = pgx_by_pop.get(pop)
        transplant_row = latest.get(iso3)

        point = EquityPoint(
            population_code=pop,
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else iso3,
            transplants_pmp=float(transplant_row["transplants_pmp"]) if transplant_row else None,
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pgx_delta_tacrolimus=float(pgx.get("delta_vs_baseline", 0)) if pgx else 0.0,
            pct_requiring_change=float(pgx.get("percentage_requiring_change", 0)) if pgx else 0.0,
        )
        point.compute(global_median)
        points.append(point)

    return sorted(points, key=lambda p: p.equity_score, reverse=True)
