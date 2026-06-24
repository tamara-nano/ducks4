#!/usr/bin/env python3

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Optional

import pandas as pd


LOG = logging.getLogger("curate_annotation_bed")


# repeat size thresholds
FULL_MIN = 3000
FULL_MAX = 3600
S_MIN    = 1500
S_MAX    = 2999
INTERNAL_MIN_BP = 1500

# structural gap thresholds
GAP_LARGE_THRESHOLD = 10000   # gap >10kb between D4Z4 hits = structural gap
GAP_SMALL_MAX       = 5000    # gap between real neighbouring RUs should be <5kb
CTRL_MAX_LEN        = 2000    # pre-gap chr4_ctrl fragments are short <2kb

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    )


# BED / bedGraph 
def read_bed(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", header=None, comment="#")
    if df.shape[1] < 4:
        raise ValueError(f"{path} does not look like a 4-column BED.")

    df = df.iloc[:, :4].copy()
    df.columns = ["contig", "start", "end", "raw_name"]
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"]   = pd.to_numeric(df["end"],   errors="coerce")

    if df[["start", "end"]].isna().any().any():
        raise ValueError(f"{path} contains invalid coordinates.")

    df = df.reset_index(drop=True)
    df["row_id"]    = df.index
    df["length_bp"] = df["end"] - df["start"]
    return df


def read_bedgraph(path: Path, value_col: int = 3) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) <= value_col:
                LOG.debug("Skipping non-bedGraph line in %s: %s", path, line)
                continue
            rows.append(fields)

    if not rows:
        raise ValueError(f"{path} contains no usable bedGraph rows after filtering headers.")

    df  = pd.DataFrame(rows)
    out = df.iloc[:, [0, 1, 2, value_col]].copy()
    out.columns = ["contig", "start", "end", "bedgraph_methylation"]
    out["start"] = pd.to_numeric(out["start"], errors="coerce")
    out["end"]   = pd.to_numeric(out["end"],   errors="coerce")
    out["bedgraph_methylation"] = pd.to_numeric(out["bedgraph_methylation"], errors="coerce")

    if out[["start", "end"]].isna().any().any():
        bad = out[out["start"].isna() | out["end"].isna()]
        raise ValueError(
            f"{path} contains invalid coordinates after header filtering. "
            f"Example bad rows:\n{bad.head().to_string(index=False)}"
        )

    out = out.reset_index(drop=True)
    out["row_id"] = out.index
    return out


def validate_bed_and_bedgraph_are_mirrored(
    bed_df: pd.DataFrame,
    bedgraph_df: pd.DataFrame,
) -> bool:
    if len(bed_df) != len(bedgraph_df):
        LOG.warning(
            "BED and bedGraph row counts differ: %d BED rows vs %d bedGraph rows.",
            len(bed_df), len(bedgraph_df),
        )
        return False

    coords_match = (
        bed_df[["contig", "start", "end"]].reset_index(drop=True)
        .equals(bedgraph_df[["contig", "start", "end"]].reset_index(drop=True))
    )
    if not coords_match:
        LOG.warning("BED and bedGraph coordinates are not identical in order.")
        return False

    LOG.info("Validated mirrored BED/bedGraph rows: %d entries", len(bed_df))
    return True


def align_bedgraph_to_bed(
    bed_df: pd.DataFrame,
    bedgraph_df: pd.DataFrame,
) -> pd.DataFrame:
    left  = bed_df[["row_id", "contig", "start", "end"]].copy()
    right = bedgraph_df[["contig", "start", "end", "bedgraph_methylation"]].copy()

    merged = left.merge(right, on=["contig", "start", "end"], how="left", validate="one_to_one")

    n_missing = merged["bedgraph_methylation"].isna().sum()
    if n_missing > 0:
        missing = merged.loc[merged["bedgraph_methylation"].isna(), ["contig", "start", "end"]]
        raise ValueError(
            f"Could not align bedGraph to all BED rows. Missing {n_missing} rows. "
            f"Example:\n{missing.head().to_string(index=False)}"
        )
    return merged.reset_index(drop=True)


# structural gap analysis
def find_structural_gaps(d4z4_df: pd.DataFrame, threshold: int = GAP_LARGE_THRESHOLD) -> list[dict]:
    """
    Find large structural gaps (>threshold bp) between consecutive D4Z4_like hits.

    Returns list of dicts:
      {gap_after_idx, gap_size, pre_gap_rows, post_gap_rows}

    gap_after_idx = positional index in d4z4_df (0-based) after which the gap occurs.
    pre_gap_rows  = d4z4_df rows before gap
    post_gap_rows = d4z4_df rows after gap
    """
    if d4z4_df.empty or len(d4z4_df) < 2:
        return []

    sorted_df = d4z4_df.sort_values("start").reset_index(drop=True)
    gaps = []

    for i in range(len(sorted_df) - 1):
        gap = int(sorted_df.at[i + 1, "start"]) - int(sorted_df.at[i, "end"])
        if gap > threshold:
            gaps.append({
                "gap_after_idx": i,
                "gap_size":      gap,
                "pre_gap_rows":  sorted_df.iloc[:i + 1].copy(),
                "post_gap_rows": sorted_df.iloc[i + 1:].copy(),
            })

    return gaps


def detect_deleted_d4f104s1(df: pd.DataFrame) -> bool:
    """
    Detect alleles where D4F104S1 is deleted based on structural gap pattern.

    Expected pattern for chr4 deleted-D4F104S1 allele:
      [4qB_pLAM_lowpid]                 ← pseudo_pLAM proximal marker
      [c10 ~1550bp INVERTED]  ──┐        ← chr4_ctrl pair — inverted D4Z4
      [tiny c4 fragment ~111bp] ┘        ← chr4_ctrl overlap fragment
      [large gap >10kb]                  ← deleted D4F104S1 region
      [optional short pre-RU fragment]   ← partial overlap hit
      [full c4 RUs ≥3kb, no big gaps]
      [4qA_pLAM]
      [DUX4_end]

    The inverted c10 hit is a positive chr4 marker — chr10 would not show this.
    """
    # must have pseudo_pLAM and NO real D4F104S1
    if df[df["feature_class"] == "pseudo_pLAM"].empty:
        return False
    if not df[df["feature_class"] == "D4F104S1"].empty:
        return False

    # must have real pLAM
    true_plam = df[
        (df["feature_class"] == "pLAM") &
        ~df["raw_name"].astype(str).str.contains("lowpid|low-pid", case=False, na=False)
    ]
    if true_plam.empty:
        return False

    d4z4 = df[df["feature_class"] == "D4Z4_like"].sort_values("start").reset_index(drop=True)
    if len(d4z4) < 3:
        return False

    gaps = find_structural_gaps(d4z4, threshold=GAP_LARGE_THRESHOLD)
    if not gaps:
        return False

    # use first large gap
    gap_info  = gaps[0]
    pre_gap   = gap_info["pre_gap_rows"]
    post_gap  = gap_info["post_gap_rows"].reset_index(drop=True)

    if post_gap.empty:
        return False

    # pre-gap hits must all be short (<2kb) — chr4_ctrl inverted pair
    if not (pre_gap["length_bp"] < CTRL_MAX_LEN).all():
        return False

    # positive chr4 marker: inverted c10 hit ~1000-2000bp in pre-gap
    has_inverted_c10 = (
        pre_gap["repeat_origin"].eq("c10") &
        pre_gap["length_bp"].between(500, 2000)
    ).any()

    if not has_inverted_c10:
        LOG.debug("No chr4_ctrl inverted c10 hit — not a deleted D4F104S1 pattern")
        return False

    # post-gap: find first full RU ≥3kb
    full_ru_mask = post_gap["length_bp"] >= FULL_MIN
    if not full_ru_mask.any():
        return False

    first_full_idx = int(full_ru_mask.idxmax())
    real_rus = post_gap.iloc[first_full_idx:].reset_index(drop=True)

    if len(real_rus) < 2:
        return False

    # no large internal gaps between real RUs
    internal_gaps = find_structural_gaps(real_rus, threshold=GAP_SMALL_MAX)
    if internal_gaps:
        LOG.debug(
            "Large internal gap between post-gap RUs (%dbp) — not deleted D4F104S1 pattern",
            internal_gaps[0]["gap_size"],
        )
        return False

    LOG.info(
        "Detected deleted D4F104S1 (chr4_ctrl inverted pair confirmed): "
        "gap of %dbp, %d pre-gap fragments, %d real RUs post-gap.",
        gap_info["gap_size"],
        len(pre_gap),
        len(real_rus),
    )
    return True


# feature classification
def classify_raw_feature(raw_name: str) -> str:
    r = str(raw_name)

    if "chr10_ctrl" in r:
        return "chr10_ctrl"
    if "CLUHP4" in r:
        return "CLUHP4"
    if "DUX4_end" in r:
        return "DUX4_end"
    if "D4F104S1" in r:
        return "D4F104S1"
    # ← lowpid MUST come before general pLAM check
    if "pLAM_lowpid" in r or "PLAM_lowpid" in r:
        return "pseudo_pLAM"
    if "pLAM" in r or "PLAM" in r:
        return "pLAM"
    if "D4Z4" in r:
        return "D4Z4_like"
    if re.match(r"^RU\d+_", r):
        return "D4Z4_like"
    return "other"


def infer_repeat_origin(raw_name: str) -> Optional[str]:
    r = str(raw_name)
    if "4q35_D4Z4" in r:
        return "c4"
    if "10q26_D4Z4" in r:
        return "c10"
    return None


def infer_haplotype_anchor(raw_name: str) -> Optional[str]:
    r = str(raw_name)
    for tag in ["4qA", "4qB", "10qA"]:
        if r.startswith(f"{tag}_D4F104S1") or r.startswith(f"{tag}_pLAM"):
            return tag
    return None


def add_base_columns(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["feature_class"]    = out["raw_name"].apply(classify_raw_feature)
    out["repeat_origin"]    = out["raw_name"].apply(infer_repeat_origin)
    out["anchor_haplotype"] = out["raw_name"].apply(infer_haplotype_anchor)
    return out


# anchor detection
def find_anchor_indices(
    df: pd.DataFrame,
) -> tuple[Optional[int], Optional[int], Optional[str], Optional[str]]:
    """
    Identify the true array anchors.

    Returns: marker_idx, plam_idx, haplotype_tag, orientation

    Orientation values:
      "forward"               — D4F104S1 ... pLAM  (normal)
      "inverted"              — pLAM ... D4F104S1  (inverted read)
      "forward_no_D4F104S1"  — D4F104S1 deleted; proximal boundary = first real RU
    """
    marker_df = df[df["feature_class"] == "D4F104S1"].copy()
    plam_df   = df[
        (df["feature_class"] == "pLAM") &
        ~df["raw_name"].astype(str).str.contains("lowpid|low-pid", case=False, na=False)
    ].copy()

    # normal forward / inverted read
    if not marker_df.empty and not plam_df.empty:
        # prefer forward
        for m_idx, m_row in marker_df.iterrows():
            downstream_plam = plam_df[plam_df["start"] >= m_row["end"]]
            if not downstream_plam.empty:
                p_idx = downstream_plam.index[0]
                return int(m_idx), int(p_idx), m_row["anchor_haplotype"], "forward"

        # inverted
        for p_idx, p_row in plam_df.iterrows():
            downstream_marker = marker_df[marker_df["start"] >= p_row["end"]]
            if not downstream_marker.empty:
                m_idx = (
                    downstream_marker.index[-1]
                    if len(downstream_marker) > 1
                    else downstream_marker.index[0]
                )
                return int(m_idx), int(p_idx), downstream_marker.iloc[0]["anchor_haplotype"], "inverted"

        # fallback: positional
        m_idx = int(marker_df.index[0])
        p_idx = int(plam_df.index[0])
        hap   = marker_df.iloc[0]["anchor_haplotype"]
        orientation = "forward" if marker_df.iloc[0]["start"] < plam_df.iloc[0]["start"] else "inverted"
        return m_idx, p_idx, hap, orientation

    # deleted D4F104S1 fallback
    if detect_deleted_d4f104s1(df):
        true_plam = plam_df if not plam_df.empty else df[df["feature_class"] == "pLAM"]
        if true_plam.empty:
            return None, None, None, None

        d4z4 = df[df["feature_class"] == "D4Z4_like"].sort_values("start").reset_index(drop=True)

        gaps = find_structural_gaps(d4z4, threshold=GAP_LARGE_THRESHOLD)
        if not gaps:
            return None, None, None, None

        post_gap = gaps[0]["post_gap_rows"].reset_index(drop=True)

        # first full RU ≥3kb after the gap = proximal array boundary
        full_mask = post_gap["length_bp"] >= FULL_MIN
        if not full_mask.any():
            return None, None, None, None

        first_real_start = int(post_gap.loc[full_mask.idxmax(), "start"])
        m_candidates = df[df["start"] == first_real_start]
        if m_candidates.empty:
            return None, None, None, None

        m_idx = int(m_candidates.index[0])
        p_idx = int(true_plam.sort_values("start").index[-1])
        hap   = true_plam.at[p_idx, "anchor_haplotype"]

        return m_idx, p_idx, hap, "forward_no_D4F104S1"

    return None, None, None, None


# overlap / collapse
def overlap_fraction(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    ov = max(0, min(a_end, b_end) - max(a_start, b_start))
    if ov == 0:
        return 0.0
    shorter = min(a_end - a_start, b_end - b_start)
    return ov / shorter if shorter > 0 else 0.0


def collapse_duplicate_d4z4(
    subdf: pd.DataFrame,
    preferred_origin: Optional[str] = None,
) -> pd.DataFrame:
    """
    Collapse strongly overlapping D4Z4-like hits inside the candidate array block.
    Best BLAST hit (highest pident → bitscore → length) wins.
    """
    if subdf.empty:
        return subdf.copy()

    rows = subdf.sort_values(["start", "end"]).copy()
    rows["keep_overlap"] = True
    rows["drop_reason"]  = ""

    idxs = list(rows.index)
    used = set()

    for i_pos, i in enumerate(idxs):
        if i in used:
            continue

        ri    = rows.loc[i]
        group = [i]
        used.add(i)

        for j in idxs[i_pos + 1:]:
            if j in used:
                continue
            rj   = rows.loc[j]
            frac = overlap_fraction(int(ri["start"]), int(ri["end"]), int(rj["start"]), int(rj["end"]))
            if frac >= 0.8:
                group.append(j)
                used.add(j)

        if len(group) == 1:
            continue

        grp        = rows.loc[group].copy()
        sort_cols  = []
        ascending  = []

        if "pident" in grp.columns:
            sort_cols.append("pident");   ascending.append(False)
        if "bitscore" in grp.columns:
            sort_cols.append("bitscore"); ascending.append(False)
        sort_cols.append("length_bp");    ascending.append(False)

        grp      = grp.sort_values(by=sort_cols, ascending=ascending)
        keep_idx = grp.index[0]

        for idx in group:
            if idx != keep_idx:
                rows.at[idx, "keep_overlap"] = False
                rows.at[idx, "drop_reason"]  = "overlapping_duplicate_d4z4_lower_pid"

    return rows


# terminal repeat classification
def classify_terminal_repeat(length_bp: int, origin: str) -> str:
    if FULL_MIN <= length_bp <= FULL_MAX:
        return f"{origin}-L"
    if S_MIN <= length_bp <= S_MAX:
        return f"{origin}-S"
    if length_bp < S_MIN:
        return f"{origin}-fragment"
    return f"{origin}-L"


def choose_true_terminal_repeat(
    candidate: pd.DataFrame,
) -> tuple[int, list[int], list[int]]:
    """
    Returns: terminal_idx, true_array_indices, trailing_fragment_indices
    """
    if candidate.empty:
        raise ValueError("choose_true_terminal_repeat() received empty candidate table.")

    cand = candidate.copy()
    countable_mask = (
        cand["length_bp"].between(FULL_MIN, FULL_MAX, inclusive="both") |
        cand["length_bp"].between(S_MIN,    S_MAX,    inclusive="both")
    )

    if countable_mask.any():
        terminal_idx = cand.loc[countable_mask].index[-1]
        terminal_pos = cand.index.get_loc(terminal_idx)
    else:
        terminal_pos = len(cand) - 1
        terminal_idx = cand.index[terminal_pos]

    return (
        terminal_idx,
        cand.index[:terminal_pos + 1].tolist(),
        cand.index[terminal_pos + 1:].tolist(),
    )


# main curation

PSEUDO_PLAM_MAX_LEN = 200  # real pLAM ≥200bp; proxy pseudo-hits are shorter

def reclassify_proximal_pseudo_plam(df: pd.DataFrame) -> pd.DataFrame:

    out      = df.copy()
    plam_hits = out[out["feature_class"] == "pLAM"].copy()

    if len(plam_hits) < 2:
        return out

    short_plam = plam_hits[plam_hits["length_bp"] < PSEUDO_PLAM_MAX_LEN]
    long_plam  = plam_hits[plam_hits["length_bp"] >= PSEUDO_PLAM_MAX_LEN]

    if short_plam.empty or long_plam.empty:
        return out

    for idx in short_plam.index:
        out.at[idx, "feature_class"]   = "pseudo_pLAM"
        out.at[idx, "curated_type"]    = "pseudo_pLAM"
        out.at[idx, "curation_reason"] = "short_plam_reclassified_as_pseudo_lowpid"
        LOG.info(
            "Reclassified short pLAM hit at %d-%d (%dbp) as pseudo_pLAM "
            "(likely 4qB_pLAM_lowpid missing suffix from parse_blast_to_bed).",
            int(out.at[idx, "start"]),
            int(out.at[idx, "end"]),
            int(out.at[idx, "length_bp"]),
        )

    return out
  
def curate_annotation(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy().sort_values(["start", "end"]).reset_index(drop=True)

    out["curated_name"]    = out["raw_name"]
    out["curated_type"]    = out["feature_class"]
    out["is_array_member"] = False
    out["count_for_RU"]    = False
    out["ru_index"]        = pd.NA
    out["is_distal"]       = False
    out["curation_reason"] = ""

    # reclassify short pLAM hits missed by parse_blast_to_bed 
    out = reclassify_proximal_pseudo_plam(out)

    # pseudo_pLAM: keep visible
    pseudo_plam_mask = out["feature_class"] == "pseudo_pLAM"
    out.loc[pseudo_plam_mask, "curated_name"]    = out.loc[pseudo_plam_mask, "raw_name"]
    out.loc[pseudo_plam_mask, "curated_type"]    = "pseudo_pLAM"
    out.loc[pseudo_plam_mask, "curation_reason"] = "low_identity_plam_pseudohit"


    marker_idx, plam_idx, hap, orientation = find_anchor_indices(out)

    if marker_idx is None or plam_idx is None:
        LOG.warning("Could not identify anchor pair. Returning lightly classified table only.")
        for fc, cn, ct, cr in [
            ("CLUHP4",     "CLUHP4",     "CLUHP4",     "framing_control"),
            ("DUX4_end",   "DUX4_end",   "DUX4_end",   "framing_control"),
            ("chr10_ctrl", "chr10_ctrl", "chr10_ctrl",  "chromosome10_control"),
        ]:
            mask = out["feature_class"] == fc
            out.loc[mask, "curated_name"]    = cn
            out.loc[mask, "curated_type"]    = ct
            out.loc[mask, "curation_reason"] = cr
        out.loc[out["curation_reason"] == "", "curation_reason"] = "unmodified"
        return out

    expected_origin = None
    if hap in {"4qA", "4qB"}:
        expected_origin = "c4"
    elif hap == "10qA":
        expected_origin = "c10"

    # framing controls
    for fc, cn, ct in [
        ("CLUHP4",     "CLUHP4",     "CLUHP4"),
        ("DUX4_end",   "DUX4_end",   "DUX4_end"),
        ("chr10_ctrl", "chr10_ctrl", "chr10_ctrl"),
    ]:
        out.loc[out["feature_class"] == fc, "curated_name"] = cn
        out.loc[out["feature_class"] == fc, "curated_type"] = ct

    # handle deleted D4F104S1 - DBED
    effective_orientation = orientation

    if orientation == "forward_no_D4F104S1":
        effective_orientation = "forward"

        # this is the inverted D4Z4 unit on chr4
        pre_gap_hits = out[
            (out["feature_class"] == "D4Z4_like") &
            (out.index < marker_idx)
        ].copy()

        if not pre_gap_hits.empty:
            # longest pre-gap hit = the inverted D4Z4 unit → chr4_ctrl
            inverted_hit_idx = pre_gap_hits["length_bp"].idxmax()
            out.at[inverted_hit_idx, "curated_name"]    = "chr4_ctrl"
            out.at[inverted_hit_idx, "curated_type"]    = "chr4_ctrl"
            out.at[inverted_hit_idx, "curation_reason"] = (
                "chr4_ctrl_inverted_D4Z4_proximal_to_deleted_D4F104S1"
            )

            # remaining pre-gap hits: check proximity to first real RU
            # hits immediately adjacent to marker_idx = deletion boundary fragment
            # hits far from marker_idx = chr4_ctrl (fall through to outside_mask)
            first_ru_start = int(out.at[marker_idx, "start"])

            for idx in pre_gap_hits.index:
                if idx == inverted_hit_idx:
                    continue  # already labeled

                hit_end = int(out.at[idx, "end"])
                gap_to_ru1 = first_ru_start - hit_end

                if gap_to_ru1 <= 1000:
                    # immediately adjacent to RU01 — deletion boundary fragment
                    out.at[idx, "curated_name"]    = "D4Z4_fragment"
                    out.at[idx, "curated_type"]    = "internal_fragment"
                    out.at[idx, "curation_reason"] = (
                        "D4Z4_fragment_at_D4F104S1_deletion_boundary"
                    )
                # else: fall through to outside_mask → chr4_ctrl normally

        out.at[marker_idx, "curation_reason"] = "array_proximal_boundary_D4F104S1_deleted"

        LOG.info(
            "Curation: deleted D4F104S1 allele — "
            "first real RU at %d used as proximal boundary.",
            int(out.at[marker_idx, "start"]),
        )

    else:
        # normal anchor labeling
        out.at[marker_idx, "curated_name"]    = out.at[marker_idx, "raw_name"]
        out.at[marker_idx, "curated_type"]    = "D4F104S1"
        out.at[marker_idx, "curation_reason"] = "array_anchor"

    # pLAM anchor
    out.at[plam_idx, "curated_name"]    = out.at[plam_idx, "raw_name"]
    out.at[plam_idx, "curated_type"]    = "pLAM"
    out.at[plam_idx, "is_distal"]       = True
    out.at[plam_idx, "curation_reason"] = "array_terminal_anchor"

    # candidate repeat array
    left_idx  = min(marker_idx, plam_idx)
    right_idx = max(marker_idx, plam_idx)

    candidate_mask = (
        (out.index > left_idx) &
        (out.index < right_idx) &
        (out["feature_class"] == "D4Z4_like")
    )

    candidate = out.loc[candidate_mask].copy()
    candidate = collapse_duplicate_d4z4(candidate, preferred_origin=expected_origin)

    if "keep_overlap" in candidate.columns:
        candidate = candidate[candidate["keep_overlap"] == True].copy()

    kept_candidate_idx    = set(candidate.index.tolist())
    dropped_candidate_idx = set(out.loc[candidate_mask].index.tolist()) - kept_candidate_idx

    if dropped_candidate_idx:
        out.loc[list(dropped_candidate_idx), "curation_reason"] = "overlapping_duplicate_d4z4_dropped"

    # D4Z4_like hits outside the true array = chr4_ctrl
    outside_mask = (out["feature_class"] == "D4Z4_like") & ~candidate_mask
    # but don't relabel pre-gap hits already labeled in deleted-D4F case
    if orientation == "forward_no_D4F104S1":
        outside_mask = outside_mask & (out["curation_reason"] == "")
    out.loc[outside_mask, "curated_name"]    = "chr4_ctrl"
    out.loc[outside_mask, "curated_type"]    = "chr4_ctrl"
    out.loc[outside_mask, "curation_reason"] = "pseudo_d4z4_outside_true_array"

    # assign RU numbers
    if not candidate.empty:
        if effective_orientation == "forward":
            candidate = candidate.sort_values(["start", "end"]).copy()
        else:
            candidate = candidate.sort_values(["start", "end"], ascending=[False, False]).copy()

        terminal_idx, true_array_indices, trailing_fragment_indices = (
            choose_true_terminal_repeat(candidate)
        )

        # internal RUs
        for idx in true_array_indices[:-1]:
            origin    = candidate.at[idx, "repeat_origin"] or expected_origin or "c4"
            length_bp = int(candidate.at[idx, "length_bp"])

            if length_bp >= INTERNAL_MIN_BP:
                out.at[idx, "curated_name"]    = origin
                out.at[idx, "curated_type"]    = "D4Z4"
                out.at[idx, "is_array_member"] = True
                out.at[idx, "count_for_RU"]    = True
                out.at[idx, "curation_reason"] = (
                    "internal_array_repeat"
                    if effective_orientation == "forward"
                    else "internal_array_repeat_inverted"
                )
            else:
                out.at[idx, "curated_name"]    = "internal_fragment"
                out.at[idx, "curated_type"]    = "internal_fragment"
                out.at[idx, "is_array_member"] = False
                out.at[idx, "count_for_RU"]    = False
                out.at[idx, "is_distal"]       = False
                out.at[idx, "curation_reason"] = "short_internal_d4z4_fragment"

        # terminal RU
        terminal_origin = candidate.at[terminal_idx, "repeat_origin"] or expected_origin or "c4"
        terminal_len    = int(candidate.at[terminal_idx, "length_bp"])
        distal_label    = classify_terminal_repeat(terminal_len, terminal_origin)

        out.at[terminal_idx, "curated_name"]    = distal_label
        out.at[terminal_idx, "curated_type"]    = "D4Z4-L" if distal_label.endswith("-L") else "D4Z4-S"
        out.at[terminal_idx, "is_array_member"] = True
        out.at[terminal_idx, "count_for_RU"]    = True
        out.at[terminal_idx, "is_distal"]       = True
        out.at[terminal_idx, "curation_reason"] = (
            "terminal_array_repeat"
            if effective_orientation == "forward"
            else "terminal_array_repeat_inverted"
        )

        # RU numbering
        countable_ru_indices = [
            idx for idx in true_array_indices[:-1]
            if int(candidate.at[idx, "length_bp"]) >= INTERNAL_MIN_BP
        ]
        countable_ru_indices.append(terminal_idx)

        for n, idx in enumerate(countable_ru_indices, start=1):
            out.at[idx, "ru_index"] = n

        # trailing fragments
        for idx in trailing_fragment_indices:
            out.at[idx, "curated_name"]    = "distal_fragment"
            out.at[idx, "curated_type"]    = "distal_fragment"
            out.at[idx, "is_array_member"] = False
            out.at[idx, "count_for_RU"]    = False
            out.at[idx, "is_distal"]       = False
            out.at[idx, "curation_reason"] = "post_terminal_d4z4_fragment"

    # fill empty curation reasons
    for fc, cr in [
        ("CLUHP4",     "framing_control"),
        ("DUX4_end",   "framing_control"),
        ("chr10_ctrl", "chromosome10_control"),
    ]:
        mask = (out["feature_class"] == fc) & (out["curation_reason"] == "")
        out.loc[mask, "curation_reason"] = cr

    out.loc[out["curation_reason"] == "", "curation_reason"] = "unmodified"
    return out


# output
def make_display_name(row: pd.Series) -> str:
    curated_name = str(row["curated_name"])
    curated_type = str(row.get("curated_type", ""))

    if curated_type == "pseudo_pLAM":
        return f"{curated_name}_lowpid"   # ← makes it explicit in BED

    if pd.notna(row.get("ru_index")):
        ru = int(row["ru_index"])
        return f"RU{ru:02d}_{curated_name}"
    return curated_name


def make_curated_bed(curated: pd.DataFrame) -> pd.DataFrame:
    keep_types = {
        "chr4_ctrl", "chr10_ctrl",
        "CLUHP4", "D4F104S1",
        "D4Z4", "D4Z4-L", "D4Z4-S",
        "internal_fragment", "distal_fragment",
        "pLAM", "pseudo_pLAM", "DUX4_end",
    }

    bed = curated[curated["curated_type"].isin(keep_types)].copy()
    bed = bed.sort_values(["start", "end"]).reset_index(drop=True)
    bed["display_name"] = bed.apply(make_display_name, axis=1)

    bed_out = bed[[
        "row_id", "contig", "start", "end",
        "display_name", "curated_name", "curated_type",
        "ru_index", "count_for_RU", "is_array_member", "is_distal",
        "curation_reason",
    ]].copy()

    # distal_unit aggregate region
    count_for_ru = (
        curated["count_for_RU"].astype(str).str.strip().str.upper().isin(["TRUE", "WAHR", "1"])
        if "count_for_RU" in curated.columns
        else pd.Series(False, index=curated.index)
    )
    is_distal = (
        curated["is_distal"].astype(str).str.strip().str.upper().isin(["TRUE", "WAHR", "1"])
        if "is_distal" in curated.columns
        else pd.Series(False, index=curated.index)
    )

    distal_mask = (count_for_ru & is_distal) | curated["curated_type"].eq("pLAM")
    distal_sub  = curated[distal_mask].copy()

    if not distal_sub.empty:
        distal_row = {col: pd.NA for col in bed_out.columns}
        distal_row.update({
            "row_id":          -1,
            "contig":          str(distal_sub["contig"].iloc[0]),
            "start":           int(distal_sub["start"].min()),
            "end":             int(distal_sub["end"].max()),
            "display_name":    "distal_unit",
            "curated_name":    "distal_unit",
            "curated_type":    "distal_unit",
            "ru_index":        pd.NA,
            "count_for_RU":    False,
            "is_array_member": False,
            "is_distal":       True,
            "curation_reason": "aggregate_distal_repeat_plus_pLAM",
        })
        bed_out = pd.concat([bed_out, pd.DataFrame([distal_row])], ignore_index=True)
        bed_out = bed_out.sort_values(["start", "end", "display_name"]).reset_index(drop=True)

    return bed_out


def make_curated_bedgraph(
    curated_bed: pd.DataFrame,
    raw_bedgraph: pd.DataFrame,
) -> pd.DataFrame:
    keep_ids = set(curated_bed["row_id"].tolist())
    bg       = raw_bedgraph[raw_bedgraph["row_id"].isin(keep_ids)].copy()
    meta     = curated_bed[["row_id", "contig", "start", "end"]].copy()

    bg = meta.merge(
        bg[["row_id", "bedgraph_methylation"]],
        on="row_id", how="left", validate="one_to_one",
    )
    bg = bg.sort_values(["start", "end"]).reset_index(drop=True)
    return bg[["contig", "start", "end", "bedgraph_methylation"]].copy()


def make_fused_output_table(
    curated: pd.DataFrame,
    curated_bed: pd.DataFrame,
    raw_bedgraph: pd.DataFrame,
) -> pd.DataFrame:
    keep_ids = set(curated_bed["row_id"].tolist())
    bg       = raw_bedgraph[raw_bedgraph["row_id"].isin(keep_ids)].copy()
    fused    = curated[curated["row_id"].isin(keep_ids)].copy()

    fused = fused.merge(
        bg[["row_id", "bedgraph_methylation"]],
        on="row_id", how="left", validate="one_to_one",
    )
    fused = fused.sort_values(["start", "end"]).reset_index(drop=True)
    fused["display_name"] = fused.apply(make_display_name, axis=1)

    cols = [
        "row_id", "contig", "start", "end",
        "raw_name", "display_name", "length_bp",
        "feature_class", "repeat_origin", "anchor_haplotype",
        "curated_name", "curated_type",
        "ru_index", "count_for_RU", "is_array_member", "is_distal",
        "curation_reason", "bedgraph_methylation",
    ]
    return fused[cols].copy()

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Curate raw allele annotation BED into per-element anchored FSHD array annotation. "
            "Handles normal arrays, inverted reads, and deleted D4F104S1 alleles. "
            "Optionally applies the same keep/drop decisions to a matching bedGraph."
        )
    )
    p.add_argument("--bed",           required=True,  help="Input raw 4-column BED")
    p.add_argument("--bedgraph",                      help="Input mirrored annotation-level bedGraph")
    p.add_argument("--bedgraph-value-col", type=int, default=3,
                   help="0-based column index of methylation value in bedGraph (default: 3)")
    p.add_argument("--out-curation",  required=True,  help="Detailed curated TSV")
    p.add_argument("--out-bed",       required=True,  help="Curated BED for IGV")
    p.add_argument("--out-fused",                     help="Curated annotation + methylation TSV")
    p.add_argument("--out-bedgraph",                  help="Curated bedGraph for IGV")
    p.add_argument("--verbose",       action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    setup_logging(args.verbose)

    df = read_bed(Path(args.bed))
    df = add_base_columns(df)

    raw_bedgraph = None
    if args.bedgraph:
        raw_bedgraph = read_bedgraph(Path(args.bedgraph), value_col=args.bedgraph_value_col)
        mirrored_ok  = validate_bed_and_bedgraph_are_mirrored(df, raw_bedgraph)
        if not mirrored_ok:
            LOG.warning("BED and bedGraph not perfectly mirrored — falling back to coordinate alignment.")
            raw_bedgraph = align_bedgraph_to_bed(df, raw_bedgraph)

    curated     = curate_annotation(df)
    curated_bed = make_curated_bed(curated)

    Path(args.out_curation).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_bed).parent.mkdir(parents=True, exist_ok=True)

    curated.to_csv(args.out_curation, sep="\t", index=False)

    curated_bed[["contig", "start", "end", "display_name"]].to_csv(
        args.out_bed, sep="\t", index=False, header=False
    )

    LOG.info("Wrote curation table: %s", args.out_curation)
    LOG.info("Wrote curated BED:    %s", args.out_bed)

    if raw_bedgraph is not None:
        if args.out_bedgraph:
            Path(args.out_bedgraph).parent.mkdir(parents=True, exist_ok=True)
            make_curated_bedgraph(curated_bed, raw_bedgraph).to_csv(
                args.out_bedgraph, sep="\t", index=False, header=False
            )
            LOG.info("Wrote curated bedGraph: %s", args.out_bedgraph)

        if args.out_fused:
            Path(args.out_fused).parent.mkdir(parents=True, exist_ok=True)
            make_fused_output_table(curated, curated_bed, raw_bedgraph).to_csv(
                args.out_fused, sep="\t", index=False
            )
            LOG.info("Wrote fused table: %s", args.out_fused)

    elif args.out_bedgraph or args.out_fused:
        LOG.warning(
            "--out-bedgraph/--out-fused requested but no --bedgraph provided. "
            "Only annotation outputs were written."
        )


if __name__ == "__main__":
    main()
