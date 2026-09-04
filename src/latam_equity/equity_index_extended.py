"""Extended LATAM Equity Index using gnomAD v3 HWE methodology for all countries.

All PGx metrics are computed consistently via Hardy-Weinberg from gnomAD r3
allele frequencies (see ingest/gnomad_pgx.py --multi). This allows direct
comparison across all countries on the same methodological scale.

PGx gap = mean of |delta_vs_ceu| across 5 focal CPIC Level A genes:
  CYP3A5  / tacrolimus   — expressor framing: delta = (p_expr_pop - p_expr_ceu) * 100
  CYP2C19 / clopidogrel  — LOF carrier framing: delta = (p_carrier_pop - p_carrier_ceu) * 100
  VKORC1  / warfarin     — LOF carrier framing
  DPYD    / fluorouracil — LOF carrier framing
  CYP2C9  / fluvastatin  — LOF carrier framing

For CYP3A5 specifically, p_expressor = 1 - AF(*3)^2 is used (not p_lof_carrier)
because the clinically relevant gap is in expressor rate, not *3 carrier rate.
Both metrics are stored on the EquityPointExtended object for transparency.

Population → Country mapping:
  1000G pops (gnomAD 1kg:*): MXL→MEX, PEL→PER, CLM→COL, PUR→PRI
  gnomAD AMR (proxy):         AMR→ARG, BRA, CHL, URY, BOL, CRI, ECU,
                                       GTM, HND, NIC, PAN, PRY, SLV, VEN
"""
from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from src.latam_equity.equity_index import ISO3_TO_COUNTRY_NAME, PRI_TRANSPLANTS_PMP_ESTIMATE

# 1000G populations with direct country mapping
LATAM_1KG: dict[str, str] = {
    "1kg:mxl": "MEX",
    "1kg:pel": "PER",
    "1kg:clm": "COL",
    "1kg:pur": "PRI",
}

# Countries using gnomAD AMR as PGx proxy (no direct 1000G population match)
GNOMAD_AMR_COUNTRIES: list[str] = [
    "ARG", "BRA", "CHL", "URY",   # Southern Cone
    "BOL", "ECU", "PRY", "VEN",   # Andean / Gran Chaco
    "CRI", "GTM", "HND", "NIC",   # Central America
    "PAN", "SLV",                  # Central America cont.
]

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
    # CYP3A5-specific (expressor framing, for backward compat)
    pct_expressor: float           # % with ≥1 functional CYP3A5*1 allele (gnomAD HWE)
    delta_vs_ceu_pp: float         # CYP3A5 expressor delta from CEU baseline in pp
    # Composite PGx gap across 5 CPIC genes
    pgx_composite_gap: float       # mean |delta| across all focal genes
    pgx_per_gene: dict[str, float] # {gene: abs_delta_pp} breakdown
    pgx_data_source: str           # "gnomad_r3_hwe" or "gnomad_r3_hwe_amr_proxy"
    pgx_is_proxy: bool             # True = AMR proxy, not direct population match
    transplants_pmp_is_estimate: bool = False
    access_gap: float = field(init=False)
    pgx_gap: float = field(init=False)      # = pgx_composite_gap
    equity_score: float = field(init=False)

    def compute(self, global_median_pmp: float) -> None:
        pmp = self.transplants_pmp or 0.0
        self.access_gap = global_median_pmp - pmp
        self.pgx_gap = self.pgx_composite_gap
        self.equity_score = round(self.access_gap + self.pgx_gap, 3)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _index_by_pop(gnomad_multi: dict) -> dict[str, dict[str, dict]]:
    """
    Build index: {gene: {gnomad_pop_id: row_dict}}.
    Input is gnomad_pgx_multi.json (fetch_all_focal_variants output).
    """
    result: dict[str, dict[str, dict]] = {}
    for gene, gene_data in gnomad_multi.get("genes", {}).items():
        result[gene] = {row["gnomad_pop_id"]: row for row in gene_data.get("populations", [])}
    return result


def _expressor_delta(gene_pops: dict[str, dict], pop_id: str, ceu_pop_id: str = "1kg:ceu") -> tuple[float, float]:
    """
    For CYP3A5: compute expressor-framing delta from af_lof stored in multi-gene data.
    Returns (pct_expressor, delta_vs_ceu_pp).
    """
    ceu_row = gene_pops.get(ceu_pop_id, {})
    af_ceu = ceu_row.get("af_lof", 0.0)
    p_expr_ceu = 1.0 - af_ceu ** 2

    pop_row = gene_pops.get(pop_id, {})
    af_pop = pop_row.get("af_lof", af_ceu)
    p_expr_pop = 1.0 - af_pop ** 2

    return round(p_expr_pop * 100, 4), round((p_expr_pop - p_expr_ceu) * 100, 4)


def _composite_pgx(
    gene_index: dict[str, dict[str, dict]],
    pop_id: str,
    focal_genes: tuple[str, ...] = ("CYP3A5", "CYP2C19", "VKORC1", "DPYD", "CYP2C9"),
) -> tuple[float, float, dict[str, float]]:
    """
    Compute composite PGx gap for a population across all focal genes.

    CYP3A5 uses expressor framing; all others use delta_vs_ceu_pp from the
    LOF-carrier computation in gnomad_pgx_multi.json.

    Returns:
        (pct_expressor_cyp3a5, delta_cyp3a5_pp, per_gene_abs_delta_dict, composite_gap)
    Actually returns (pct_expressor, delta_cyp3a5_pp, per_gene, composite_gap)
    packed as a 4-tuple.
    """
    per_gene: dict[str, float] = {}

    # CYP3A5: expressor framing
    cyp3a5_pops = gene_index.get("CYP3A5", {})
    pct_expressor, delta_cyp3a5 = _expressor_delta(cyp3a5_pops, pop_id)
    if cyp3a5_pops.get(pop_id) is not None:
        per_gene["CYP3A5"] = round(abs(delta_cyp3a5), 4)

    # Other genes: LOF-carrier framing (delta_vs_ceu_pp stored in file)
    for gene in focal_genes:
        if gene == "CYP3A5":
            continue
        pops = gene_index.get(gene, {})
        row = pops.get(pop_id)
        if row is not None and row.get("delta_vs_ceu_pp") is not None:
            per_gene[gene] = round(abs(float(row["delta_vs_ceu_pp"])), 4)

    composite = round(sum(per_gene.values()) / len(per_gene), 4) if per_gene else 0.0
    return pct_expressor, delta_cyp3a5, per_gene, composite


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_equity_points_extended(
    world_transplants: list[dict],
    gnomad_multi: dict,
) -> list[EquityPointExtended]:
    """
    Build equity points for all LATAM countries using gnomAD multi-gene HWE PGx metrics.

    Args:
        world_transplants: rows from world-transplants.json
        gnomad_multi: parsed gnomad_pgx_multi.json (from ingest/gnomad_pgx.py --multi)
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

    gene_index = _index_by_pop(gnomad_multi)
    points: list[EquityPointExtended] = []

    # 1000G-mapped countries
    for gnomad_pop_id, iso3 in LATAM_1KG.items():
        transplant_row = latest.get(iso3)
        use_pri_estimate = iso3 == "PRI" and transplant_row is None

        pct_expr, delta_cyp3a5, per_gene, composite = _composite_pgx(gene_index, gnomad_pop_id)

        pt = EquityPointExtended(
            population_code=gnomad_pop_id.replace("1kg:", "").upper(),
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else ISO3_TO_COUNTRY_NAME.get(iso3, iso3),
            transplants_pmp=(
                PRI_TRANSPLANTS_PMP_ESTIMATE if use_pri_estimate
                else float(transplant_row["transplants_pmp"]) if transplant_row else None
            ),
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pct_expressor=pct_expr,
            delta_vs_ceu_pp=delta_cyp3a5,
            pgx_composite_gap=composite,
            pgx_per_gene=per_gene,
            pgx_data_source="gnomad_r3_hwe",
            pgx_is_proxy=False,
            transplants_pmp_is_estimate=use_pri_estimate,
        )
        pt.compute(global_median)
        points.append(pt)

    # gnomAD AMR-proxy countries
    pct_expr_amr, delta_amr, per_gene_amr, composite_amr = _composite_pgx(gene_index, "amr")
    for iso3 in GNOMAD_AMR_COUNTRIES:
        transplant_row = latest.get(iso3)
        pt = EquityPointExtended(
            population_code=AMR_PROXY_POP_CODE,
            country_iso3=iso3,
            country_name=transplant_row.get("country_name", iso3) if transplant_row else ISO3_TO_COUNTRY_NAME.get(iso3, iso3),
            transplants_pmp=float(transplant_row["transplants_pmp"]) if transplant_row else None,
            total_transplants=int(transplant_row.get("total_transplants") or 0) if transplant_row else None,
            data_year=int(transplant_row.get("year") or 0) if transplant_row else None,
            pct_expressor=pct_expr_amr,
            delta_vs_ceu_pp=delta_amr,
            pgx_composite_gap=composite_amr,
            pgx_per_gene=per_gene_amr,
            pgx_data_source="gnomad_r3_hwe_amr_proxy",
            pgx_is_proxy=True,
            transplants_pmp_is_estimate=False,
        )
        pt.compute(global_median)
        points.append(pt)

    return sorted(points, key=lambda p: p.equity_score, reverse=True)
