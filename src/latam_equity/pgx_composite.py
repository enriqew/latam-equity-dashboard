"""Compute composite PGx equity gap across multiple CPIC Level A gene-drug pairs.

Rationale for gene-drug selection (transplant context):
  - tacrolimus / CYP3A5   : primary calcineurin-inhibitor immunosuppressant
  - clopidogrel / CYP2C19 : antiplatelet for post-transplant cardiovascular events
  - warfarin    / VKORC1  : anticoagulation post-transplant
  - fluorouracil/ DPYD    : chemotherapy for malignancies in immunosuppressed patients
  - fluvastatin / CYP2C9  : statin for cardiovascular protection post-transplant

Excluded:
  - TPMT: delta_vs_baseline = 0 for all LATAM populations → no population signal
  - G6PD: CEU baseline = 100% while LATAM = 1-8%, likely a pipeline artefact
           that would dominate and distort the composite
"""
from __future__ import annotations

# One clinically representative drug per gene for transplant patients.
# All must be CPIC Level A.
CPIC_FOCAL_PAIRS: list[tuple[str, str]] = [
    ("tacrolimus",   "CYP3A5"),
    ("clopidogrel",  "CYP2C19"),
    ("warfarin",     "VKORC1"),
    ("fluorouracil", "DPYD"),
    ("fluvastatin",  "CYP2C9"),
]


def compute_composite_gap(
    drug_impact: list[dict],
    population_code: str,
) -> tuple[float, dict[str, float]]:
    """
    Compute the composite PGx gap for a population.

    For each focal gene-drug pair, look up abs(delta_vs_baseline) for the
    given population.  Returns the mean of available values plus a
    per-gene breakdown for transparency.

    Args:
        drug_impact: rows from drug_impact_summary.json
        population_code: e.g. "MXL", "PEL"
    Returns:
        (composite_gap, per_gene_breakdown)
        composite_gap: mean abs(delta) across all pairs found, or 0.0 if none
        per_gene_breakdown: {gene: abs_delta}
    """
    index: dict[tuple[str, str], dict] = {
        (r.get("drug_name", ""), r.get("gene_symbol", "")): r
        for r in drug_impact
        if r.get("population_code") == population_code
    }

    per_gene: dict[str, float] = {}
    for drug, gene in CPIC_FOCAL_PAIRS:
        row = index.get((drug, gene))
        if row is not None:
            per_gene[gene] = round(abs(float(row.get("delta_vs_baseline", 0))), 4)

    if not per_gene:
        return 0.0, per_gene

    composite = round(sum(per_gene.values()) / len(per_gene), 4)
    return composite, per_gene
