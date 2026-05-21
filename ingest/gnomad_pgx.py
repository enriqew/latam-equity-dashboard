"""Fetch allele frequencies from gnomAD v3 for PGx variants and compute carrier probabilities.

CYP3A5 determines tacrolimus metabolizer status:
  - *1 allele (REF T on + strand): functional → expressor, higher clearance
  - *3 allele (ALT C on + strand, rs776746): non-functional → non-expressor

  p_expressor is computed via Hardy-Weinberg:
      P(at least one *1) = 1 - AF(*3)^2

  delta_vs_ceu_pp = (p_expressor_pop - p_expressor_ceu) * 100
    Positive → more expressors than European baseline → more patients need higher doses

All other genes (CYP2C19, VKORC1, DPYD, CYP2C9): the ALT allele is the loss-of-function /
reduced-function / high-risk allele.  Carrier probability is computed as:
    p_lof_carrier = 1 - (1 - AF_alt)^2     (≥1 LOF allele under HWE)

For the multi-gene output (gnomad_pgx_multi.json), CYP3A5 is also expressed as p_lof_carrier
so that delta_vs_ceu_pp is directionally consistent across all genes:
    delta_vs_ceu_pp > 0 → more LOF carriers than CEU baseline
    (This inverts the sign relative to the legacy gnomad_cyp3a5.json expressor metric.)

Usage:
    python -m ingest.gnomad_pgx            # backward-compat: writes gnomad_cyp3a5.json
    python -m ingest.gnomad_pgx --multi    # writes gnomad_pgx_multi.json (see bottom of file)
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from typing import Any

import requests

from ingest._common import RAW_DIR, ensure_dirs, log

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
DATASET = "gnomad_r3"

# CYP3A5*3 (rs776746) in GRCh38: ALT=C is the *3 non-functional allele
CYP3A5_STAR3_VARIANT = "7-99672916-T-C"

# Existing LATAM 1000G pops (for validation) + gnomAD AMR aggregate (new countries)
TARGET_POPULATIONS = {"1kg:ceu", "1kg:mxl", "1kg:clm", "1kg:pel", "1kg:pur", "amr"}

# ---------------------------------------------------------------------------
# FOCAL_VARIANTS: CPIC Level A PGx variants.
#
# Variant IDs in gnomAD format: chr-pos-ref-alt (GRCh38 forward strand).
# Sources:
#   CYP3A5*3  (rs776746)  : 7-99672916-T-C   — original project variant
#   CYP2C19*2 (rs4244285) : 10-94781859-G-A  — NC_000010.11:g.94781859G>A (dbSNP)
#   VKORC1    (rs9923231) : 16-31096368-C-T  — NC_000016.10:g.31096368C>T (dbSNP)
#                            VKORC1 is on the minus strand; T (fwd) = A in cDNA notation
#                            The T allele reduces VKORC1 expression → warfarin sensitivity
#   DPYD*2A   (rs3918290) : 1-97450058-C-T   — NC_000001.11:g.97450058C>T (dbSNP)
#                            Splice-donor disruption → no DPD activity → 5-FU toxicity
#   CYP2C9*3  (rs1057910) : 10-94981296-A-C  — NC_000010.11:g.94981296A>C (dbSNP)
#                            Ile359Leu, markedly reduced enzyme activity
#
# lof_is_alt=True means ALT is the loss-of-function / reduced-function allele.
# All five variants have lof_is_alt=True.
# ---------------------------------------------------------------------------

FOCAL_VARIANTS: dict[str, dict[str, Any]] = {
    "CYP3A5": {
        "rsid": "rs776746",
        "variant_id": "7-99672916-T-C",
        "drug": "tacrolimus",
        "star_allele": "*3",
        "lof_is_alt": True,
    },
    "CYP2C19": {
        "rsid": "rs4244285",
        "variant_id": "10-94781859-G-A",
        "drug": "clopidogrel",
        "star_allele": "*2",
        "lof_is_alt": True,
    },
    "VKORC1": {
        "rsid": "rs9923231",
        "variant_id": "16-31096368-C-T",
        "drug": "warfarin",
        "star_allele": "-1639A",
        "lof_is_alt": True,
    },
    "DPYD": {
        "rsid": "rs3918290",
        "variant_id": "1-97450058-C-T",
        "drug": "fluorouracil",
        "star_allele": "*2A",
        "lof_is_alt": True,
    },
    "CYP2C9": {
        "rsid": "rs1057910",
        "variant_id": "10-94981296-A-C",
        "drug": "fluvastatin",
        "star_allele": "*3",
        "lof_is_alt": True,
    },
}

# ---------------------------------------------------------------------------
# GraphQL helpers
# ---------------------------------------------------------------------------

_VARIANT_QUERY_TEMPLATE = """
{
  variant(variantId: "%s", dataset: %s) {
    variantId
    rsid
    genome { populations { id ac an } }
  }
}
"""


def _query_gnomad_variant(variant_id: str) -> dict:
    """POST a single-variant GraphQL query to gnomAD and return the variant node."""
    query = _VARIANT_QUERY_TEMPLATE % (variant_id, DATASET)
    resp = requests.post(
        GNOMAD_API,
        json={"query": query},
        headers={"Content-Type": "application/json"},
        timeout=30,
    )
    resp.raise_for_status()
    payload = resp.json()
    if "errors" in payload:
        raise RuntimeError(f"gnomAD API error for {variant_id!r}: {payload['errors']}")
    return payload["data"]["variant"]


# ---------------------------------------------------------------------------
# HWE carrier probability helpers
# ---------------------------------------------------------------------------

def _p_expressor(af_star3: float) -> float:
    """P(at least one functional *1 allele) under Hardy-Weinberg.

    CYP3A5-specific: ALT (*3) is non-functional.
    Expressor = individual who carries ≥1 REF (*1) allele.
    """
    return round(1.0 - af_star3 ** 2, 6)


def _p_lof_carrier(af_lof: float) -> float:
    """P(at least one LOF allele) under Hardy-Weinberg.

    For genes where ALT is the LOF/risk allele (CYP2C19, VKORC1, DPYD, CYP2C9,
    and also CYP3A5 when treated uniformly in the multi-gene output).
    """
    return round(1.0 - (1.0 - af_lof) ** 2, 6)


# ---------------------------------------------------------------------------
# Legacy single-gene function — KEPT INTACT for backward compatibility
# ---------------------------------------------------------------------------

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
    """Legacy wrapper — queries only CYP3A5*3 variant."""
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


def fetch_cyp3a5() -> dict:
    """Fetch and process CYP3A5*3 population data from gnomAD v3.

    Writes data/raw/gnomad_cyp3a5.json.
    Output uses the expressor framing (p_expressor / delta_vs_ceu_pp)
    for backward compatibility with existing downstream consumers.
    """
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


# ---------------------------------------------------------------------------
# Multi-gene fetch function
# ---------------------------------------------------------------------------

def _process_variant_populations(
    variant: dict,
    gene: str,
    meta: dict[str, Any],
) -> list[dict]:
    """Build per-population records for a single focal variant.

    For all variants, ALT is the LOF/risk allele (lof_is_alt=True).
    af_lof = AC / AN for the ALT allele.
    p_lof_carrier = 1 - (1 - af_lof)^2   [HWE carrier probability]
    delta_vs_ceu_pp is computed after processing all populations, so we
    return intermediate dicts and let the caller compute the delta.

    Returns a list of dicts with keys:
        gnomad_pop_id, af_lof, p_lof_carrier, pct_lof_carrier
    (delta_vs_ceu_pp is added by the caller).
    """
    pops = {
        row["id"]: row
        for row in variant["genome"]["populations"]
        if row["id"] in TARGET_POPULATIONS and row["an"] > 0
    }

    records = []
    for pop_id, row in sorted(pops.items()):
        af_lof = row["ac"] / row["an"]
        p_carrier = _p_lof_carrier(af_lof)
        records.append(
            {
                "gnomad_pop_id": pop_id,
                "af_lof": round(af_lof, 6),
                "p_lof_carrier": p_carrier,
                "pct_lof_carrier": round(p_carrier * 100, 4),
                # delta_vs_ceu_pp filled in below
            }
        )
        log.info(
            "  [%s] %s: AF_lof=%.4f  P(carrier)=%.2f%%",
            gene,
            pop_id,
            af_lof,
            p_carrier * 100,
        )

    # Compute delta relative to CEU
    ceu_record = next((r for r in records if r["gnomad_pop_id"] == "1kg:ceu"), None)
    if ceu_record is None:
        log.warning("[%s] CEU baseline not found; delta_vs_ceu_pp will be null", gene)
        for r in records:
            r["delta_vs_ceu_pp"] = None
    else:
        p_ceu = ceu_record["p_lof_carrier"]
        log.info("  [%s] CEU baseline: P(carrier)=%.2f%%", gene, p_ceu * 100)
        for r in records:
            r["delta_vs_ceu_pp"] = round((r["p_lof_carrier"] - p_ceu) * 100, 4)

    return records


def fetch_all_focal_variants() -> dict:
    """Fetch all FOCAL_VARIANTS from gnomAD v3 and build the multi-gene JSON.

    For every gene the ALT allele is the LOF/risk allele, and carrier probability is:
        p_lof_carrier = 1 - (1 - AF_alt)^2    (HWE)
        delta_vs_ceu_pp = (p_lof_carrier_pop - p_lof_carrier_ceu) * 100
            Positive  → more LOF carriers than CEU (European) baseline.

    NOTE for CYP3A5: in this output the metric is p_lof_carrier (non-expressor
    probability), so delta_vs_ceu_pp will be OPPOSITE in sign compared to the
    legacy gnomad_cyp3a5.json expressor metric.  This is intentional: it makes
    all five genes directionally consistent (positive = more clinically at-risk
    allele carriers than CEU).

    Writes data/raw/gnomad_pgx_multi.json and returns the result dict.
    """
    genes_out: dict[str, dict] = {}

    for gene, meta in FOCAL_VARIANTS.items():
        variant_id = meta["variant_id"]
        log.info("Querying gnomAD v3 for %s (%s / %s) …", gene, meta["rsid"], variant_id)
        try:
            variant = _query_gnomad_variant(variant_id)
        except Exception as exc:  # noqa: BLE001
            log.error("Failed to fetch %s (%s): %s", gene, variant_id, exc)
            raise

        pop_records = _process_variant_populations(variant, gene, meta)

        # Patch delta_vs_ceu_pp for any record whose delta is still missing
        # (handled inside _process_variant_populations, but be explicit)
        ceu_rec = next((r for r in pop_records if r["gnomad_pop_id"] == "1kg:ceu"), None)
        if ceu_rec is not None:
            p_ceu = ceu_rec["p_lof_carrier"]
            for r in pop_records:
                if r.get("delta_vs_ceu_pp") is None:
                    r["delta_vs_ceu_pp"] = round((r["p_lof_carrier"] - p_ceu) * 100, 4)

        genes_out[gene] = {
            "variant_id": variant.get("variantId", variant_id),
            "rsid": meta["rsid"],
            "drug": meta["drug"],
            "star_allele": meta["star_allele"],
            "populations": pop_records,
        }

    result = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "dataset": DATASET,
        "methodology": (
            "HWE p_lof_carrier = 1-(1-AF_lof)^2. "
            "AF_lof is the gnomAD v3 allele frequency of the ALT allele, "
            "which is the LOF / reduced-function / high-risk allele for all five genes. "
            "Positive delta = more LOF carriers than CEU baseline."
        ),
        "genes": genes_out,
    }

    ensure_dirs()
    out_path = RAW_DIR / "gnomad_pgx_multi.json"
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    log.info(
        "Saved %s (%d genes)",
        out_path,
        len(genes_out),
    )
    return result


# ---------------------------------------------------------------------------
# Entry points
# ---------------------------------------------------------------------------

def main() -> None:
    import sys
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    if "--multi" in sys.argv:
        fetch_all_focal_variants()
    else:
        fetch_cyp3a5()


if __name__ == "__main__":
    main()
