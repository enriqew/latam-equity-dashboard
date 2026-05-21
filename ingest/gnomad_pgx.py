"""Fetch CYP3A5*3 allele frequencies from gnomAD v3 and compute expressor probabilities.

CYP3A5 determines tacrolimus metabolizer status:
  - *1 allele (REF T on + strand): functional → expressor, higher clearance
  - *3 allele (ALT C on + strand, rs776746): non-functional → non-expressor

p_expressor is computed via Hardy-Weinberg:
    P(at least one *1) = 1 - AF(*3)^2

delta_vs_ceu_pp = (p_expressor_pop - p_expressor_ceu) * 100
  Positive → more expressors than European baseline → more patients need higher doses

Usage:
    python -m ingest.gnomad_pgx
"""
from __future__ import annotations

import json
import logging

import requests

from ingest._common import RAW_DIR, ensure_dirs, log

GNOMAD_API = "https://gnomad.broadinstitute.org/api"

# CYP3A5*3 (rs776746) in GRCh38: ALT=C is the *3 non-functional allele
CYP3A5_STAR3_VARIANT = "7-99672916-T-C"
DATASET = "gnomad_r3"

# Existing LATAM 1000G pops (for validation) + gnomAD AMR aggregate (new countries)
TARGET_POPULATIONS = {"1kg:ceu", "1kg:mxl", "1kg:clm", "1kg:pel", "1kg:pur", "amr"}

_QUERY = """
{
  variant(variantId: "%s", dataset: %s) {
    variantId
    rsid
    genome { populations { id ac an } }
  }
}
""" % (
    CYP3A5_STAR3_VARIANT,
    DATASET,
)


def _query_gnomad() -> dict:
    resp = requests.post(
        GNOMAD_API,
        json={"query": _QUERY},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"gnomAD API error: {payload['errors']}")
    return payload["data"]["variant"]


def _p_expressor(af_star3: float) -> float:
    """P(at least one functional *1 allele) under Hardy-Weinberg."""
    return round(1.0 - af_star3 ** 2, 6)


def fetch_cyp3a5() -> dict:
    """Fetch and process CYP3A5*3 population data from gnomAD v3."""
    log.info("Querying gnomAD v3 for %s (rs776746) …", CYP3A5_STAR3_VARIANT)
    variant = _query_gnomad()

    pops = {
        row["id"]: row
        for row in variant["genome"]["populations"]
        if row["id"] in TARGET_POPULATIONS and row["an"] > 0
    }

    ceu = pops.get("1kg:ceu")
    if not ceu:
        raise ValueError("CEU baseline (1kg:ceu) not found in gnomAD response")

    af_ceu = ceu["ac"] / ceu["an"]
    p_expr_ceu = _p_expressor(af_ceu)
    log.info("CEU baseline: AF(*3)=%.4f  P(expressor)=%.2f%%", af_ceu, p_expr_ceu * 100)

    records = []
    for pop_id, row in sorted(pops.items()):
        af3 = row["ac"] / row["an"]
        p_expr = _p_expressor(af3)
        records.append(
            {
                "gnomad_pop_id": pop_id,
                "ac_star3": row["ac"],
                "an": row["an"],
                "af_star3": round(af3, 6),
                "p_expressor": p_expr,
                "pct_expressor": round(p_expr * 100, 4),
                "delta_vs_ceu_pp": round((p_expr - p_expr_ceu) * 100, 4),
            }
        )
        log.info(
            "  %s: AF(*3)=%.4f  P(expr)=%.2f%%  delta=%.2f pp",
            pop_id,
            af3,
            p_expr * 100,
            (p_expr - p_expr_ceu) * 100,
        )

    result = {
        "variant_id": variant["variantId"],
        "rsid": variant.get("rsid"),
        "gene": "CYP3A5",
        "star_allele": "*3",
        "drug": "tacrolimus",
        "dataset": DATASET,
        "methodology": (
            "HWE: p_expressor = 1 - AF(*3)^2. "
            "Expressor = individual with ≥1 functional CYP3A5*1 allele. "
            "These patients have higher tacrolimus clearance and may require "
            "higher doses than European-calibrated standard protocols."
        ),
        "populations": records,
    }

    ensure_dirs()
    out_path = RAW_DIR / "gnomad_cyp3a5.json"
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info("Saved %s (%d populations)", out_path, len(records))
    return result


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    fetch_cyp3a5()


if __name__ == "__main__":
    main()
