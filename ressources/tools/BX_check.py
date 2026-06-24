#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
import pysam
from pathlib import Path


def fuzzy_count(seq, pattern, max_mismatches=1):
    count = 0
    plen = len(pattern)
    for i in range(len(seq) - plen + 1):
        window = seq[i:i+plen]
        mismatches = sum(a != b for a, b in zip(window, pattern))
        if mismatches <= max_mismatches:
            count += 1
    return count

def is_distal_s_unit(ru):
    return (
        str(ru.get("curated_type", "")).strip() == "D4Z4-S"
        or str(ru.get("curated_name", "")).strip() in ["c4-S", "c10-S", "D4Z4-S"]
    )

def get_ru_intervals_from_curated_tsv(curated_tsv_path):
    df = pd.read_csv(curated_tsv_path, sep="\t", dtype=str).fillna("")

    for col in ["start", "end", "ru_index"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    ru_df = df[
        df["count_for_RU"].astype(str).str.upper().isin(["TRUE", "1"])
    ].copy()

    if ru_df.empty:
        return []

    for col in ["curated_type", "is_distal"]:
        if col not in ru_df.columns:
            ru_df[col] = ""

    return ru_df[
        ["ru_index", "contig", "start", "end", "curated_name", "curated_type", "is_distal"]
    ].to_dict("records")


def classify_ru_bx(seq, qstart, qend, max_mismatches=1):
    BINI_AVRII = "CCTAGG"   # B+ site — BinI/AvrII — intact in Chr10_D4Z4
    XAPI_APOI  = "AAATTCC"  # X+ site — XapI/ApoI  — intact in Chr4_D4Z4
    XAPI_RC    = "GGAATTT"  # XapI/ApoI reverse complement

    region = seq[qstart - 1 : qend].upper()

    n_b_exact = region.count(BINI_AVRII)
    n_b_fuzzy = fuzzy_count(region, BINI_AVRII, max_mismatches)
    n_x_exact = region.count(XAPI_APOI) + region.count(XAPI_RC)
    n_x_fuzzy = (
        fuzzy_count(region, XAPI_APOI, max_mismatches) +
        fuzzy_count(region, XAPI_RC,   max_mismatches)
    )

    b_plus  = n_b_exact > 0
    b_minus = n_b_exact == 0 and n_b_fuzzy > 0
    x_plus  = n_x_exact > 0
    x_minus = n_x_exact == 0 and n_x_fuzzy > 0

    if b_plus and not x_plus:
        d4z4_type = "Chr10_D4Z4 (B+/X-)"
    elif not b_plus and x_plus:
        d4z4_type = "Chr4_D4Z4 (B-/X+)"
    elif not b_plus and not x_plus and (b_minus or x_minus):
        d4z4_type = "Hybrid_D4Z4 (B-/X-)"
    elif b_plus and x_plus:
        d4z4_type = "ambiguous (B+/X+)"
    else:
        d4z4_type = "unclassified"

    return {
        "B_status":      "B+" if b_plus else ("B-" if b_minus else "B?"),
        "X_status":      "X+" if x_plus else ("X-" if x_minus else "X?"),
        "D4Z4_type":     d4z4_type,
        "n_BinI_exact":  n_b_exact,
        "n_BinI_fuzzy":  n_b_fuzzy,
        "n_XapI_exact":  n_x_exact,
        "n_XapI_fuzzy":  n_x_fuzzy,
    }


def classify_read_ru_from_alignment(bam_path, ru_intervals, allele_id, out_dir):
    """
    For each read in the allele BAM, check B/X restriction sites
    at each RU's reference coordinates using the alignment.
    """
    if not ru_intervals:
        return pd.DataFrame()

    bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
    rows = []
    max_ru = max(int(r["ru_index"]) for r in ru_intervals if pd.notna(r["ru_index"]))

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.query_sequence:
            continue

        seq = read.query_sequence.upper()
        aligned_pairs = dict(read.get_aligned_pairs(matches_only=True))
        ref_to_query  = {ref: qry for qry, ref in aligned_pairs.items()}

        for ru in ru_intervals:
            ru_start = int(ru["start"])
            ru_end   = int(ru["end"])
            ru_idx   = int(ru["ru_index"]) if pd.notna(ru["ru_index"]) else None

            if ru_idx is None:
                continue

            query_positions = [
                ref_to_query[ref]
                for ref in range(ru_start, ru_end)
                if ref in ref_to_query
            ]

            if not query_positions:
                continue

            q_start = min(query_positions)
            q_end   = max(query_positions) + 1
            region  = seq[q_start:q_end]

            if len(region) < 10:
                continue

            bx = classify_ru_bx(
                seq    = region,
                qstart = 1,
                qend   = len(region),
            )
            
            if (
                is_distal_s_unit(ru)
                and bx["B_status"] in ["B-", "B?"]
                and bx["X_status"] in ["X-", "X?"]
            ):
                bx["D4Z4_type"] = "D4Z4-S"

            rows.append({
                "read.id":       read.query_name,
                "allele":        allele_id,
                "RU":            ru_idx,
                "curated_name":  ru["curated_name"],
                "curated_type":  ru.get("curated_type", ""),
                "is_distal":     ru.get("is_distal", ""),
                "ref_start":     ru_start,
                "ref_end":       ru_end,
                "region_length": len(region),
                **bx,
            })

    bam.close()
    return pd.DataFrame(rows)


def process_allele(bam_path, curated_tsv_path, out_dir, allele_id):
    """Run per-RU BX classification using alignment to custom reference."""

    ru_intervals = get_ru_intervals_from_curated_tsv(curated_tsv_path)
    if not ru_intervals:
        print(f"[WARN] No RU intervals found in {curated_tsv_path}")
        return None

    print(f"[INFO] {allele_id}: {len(ru_intervals)} RU intervals from curated TSV")

    df = classify_read_ru_from_alignment(bam_path, ru_intervals, allele_id, out_dir)

    if df.empty:
        print(f"[WARN] No reads overlapped RU intervals for {allele_id}")
        return None

    df.to_csv(
        out_dir / f"{allele_id}_BX_per_RU.csv",
        sep=";", index=False, encoding="utf-8-sig"
    )

    summary = (
        df.groupby("read.id")
        .apply(lambda s: pd.Series({
            "allele":    allele_id,
            "n_RU":      len(s),
            "RU_types":  " | ".join(s["D4Z4_type"].tolist()),
            "n_c10":     (s["D4Z4_type"] == "Chr10_D4Z4 (B+/X-)").sum(),
            "n_c4":      (s["D4Z4_type"] == "Chr4_D4Z4 (B-/X+)").sum(),
            "n_hybrid":  (s["D4Z4_type"] == "Hybrid_D4Z4 (B-/X-)").sum(),
            "n_ambig":   (s["D4Z4_type"] == "ambiguous (B+/X+)").sum(),
        }))
        .reset_index()
    )

    summary.to_csv(
        out_dir / f"{allele_id}_BX_summary.csv",
        sep=";", index=False, encoding="utf-8-sig"
    )

    return df  # return per-RU df for consensus building


def build_allele_bx_consensus(per_ru_df, allele_id, min_reads=3):
    """
    Build per-RU consensus BX classification across all reads for one allele.
    Majority vote per RU position - confidence = fraction agreeing with winner.
    """
    if per_ru_df.empty:
        return pd.DataFrame()

    rows = []

    for ru_num, ru_reads in per_ru_df.groupby("RU"):
        n_total    = len(ru_reads)
        counts     = ru_reads["D4Z4_type"].value_counts()
        top_type   = counts.index[0]
        curated_type_majority = (
            ru_reads["curated_type"].mode().iloc[0]
            if "curated_type" in ru_reads.columns and not ru_reads["curated_type"].empty
            else ""
        )
        
        if curated_type_majority == "D4Z4-S":
            top_type = "D4Z4-S"
        top_count  = counts.iloc[0]
        confidence = top_count / n_total

        b_plus  = (ru_reads["B_status"] == "B+").sum()
        b_minus = (ru_reads["B_status"] == "B-").sum()
        x_plus  = (ru_reads["X_status"] == "X+").sum()
        x_minus = (ru_reads["X_status"] == "X-").sum()

        b_consensus = "B+" if b_plus > b_minus else ("B-" if b_minus > b_plus else "B?")
        x_consensus = "X+" if x_plus > x_minus else ("X-" if x_minus > x_plus else "X?")

        flag = ""
        if n_total < min_reads:
            flag = f"low_coverage(n={n_total})"
        elif confidence < 0.6:
            flag = f"ambiguous_consensus({confidence:.0%})"

        rows.append({
            "allele":                allele_id,
            "RU":                    ru_num,
            "curated_name_majority": ru_reads["curated_name"].mode().iloc[0]
                                     if "curated_name" in ru_reads.columns
                                     and not ru_reads["curated_name"].empty else "",
            "n_reads":               n_total,
            "consensus_type":        top_type,
            "B_consensus":           b_consensus,
            "X_consensus":           x_consensus,
            "confidence":            round(confidence, 3),
            "n_c10":                 (ru_reads["D4Z4_type"] == "Chr10_D4Z4 (B+/X-)").sum(),
            "n_c4":                  (ru_reads["D4Z4_type"] == "Chr4_D4Z4 (B-/X+)").sum(),
            "n_hybrid":              (ru_reads["D4Z4_type"] == "Hybrid_D4Z4 (B-/X-)").sum(),
            "n_ambig":               (ru_reads["D4Z4_type"] == "ambiguous (B+/X+)").sum(),
            "n_unclassified":        (ru_reads["D4Z4_type"] == "unclassified").sum(),
            "curated_type_majority": curated_type_majority,
            "flag":                  flag,
        })

    return pd.DataFrame(rows)


def build_allele_array_summary(consensus_df, allele_id):
    """
    Summarize the full D4Z4 array structure for one allele from per-RU consensus.
    """
    if consensus_df.empty:
        return {}

    reliable = consensus_df[consensus_df["flag"] == ""].copy()

    type_map = {
        "Chr10_D4Z4 (B+/X-)": "c10",
        "Chr4_D4Z4 (B-/X+)":  "c4",
        "Hybrid_D4Z4 (B-/X-)":   "hybrid",
        "ambiguous (B+/X+)":  "ambiguous",
        "unclassified":        "?",
        "D4Z4-S": "D4Z4-S",
    }

    array_str = "-".join(
        type_map.get(row["consensus_type"], "?")
        for _, row in consensus_df.sort_values("RU").iterrows()
    )

    reliable_types = reliable["consensus_type"].map(type_map).value_counts()

    return {
        "allele":           allele_id,
        "n_RU_total":       len(consensus_df),
        "n_RU_reliable":    len(reliable),
        "n_c10":            reliable_types.get("c10", 0),
        "n_c4":             reliable_types.get("c4", 0),
        "n_hybrid":         reliable_types.get("hybrid", 0),
        "n_ambiguous":      reliable_types.get("ambiguous", 0),
        "array_structure":  array_str,
        "low_coverage_RUs": int((consensus_df["flag"].str.startswith("low_coverage")).sum()),
        "ambiguous_RUs":    int((consensus_df["flag"].str.startswith("ambiguous_consensus")).sum()),
    }

def find_restriction_sites_in_ru(bam_path, ru_intervals, prefix_name):
    """
    For each RU, find exact reference positions of BinI (CCTAGG) and XapI (AAATTCC) sites
    by scanning reads aligned to the custom reference.
    Returns list of dicts with site positions on the reference.
    """
    BINI   = "CCTAGG"
    XAPI   = "AAATTCC"
    XAPI_RC = "GGAATTT"

    site_rows = []
    bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if not read.query_sequence:
            continue

        seq = read.query_sequence.upper()
        aligned_pairs = read.get_aligned_pairs(matches_only=True)
        query_to_ref  = {qry: ref for qry, ref in aligned_pairs}

        for ru in ru_intervals:
            ru_start = int(ru["start"])
            ru_end   = int(ru["end"])
            ru_idx   = int(ru["ru_index"]) if pd.notna(ru["ru_index"]) else None
            if ru_idx is None:
                continue

            # get query positions overlapping this RU
            ref_to_query = {ref: qry for qry, ref in aligned_pairs if ru_start <= ref < ru_end}
            if not ref_to_query:
                continue

            q_start = min(ref_to_query.values())
            q_end   = max(ref_to_query.values()) + 1
            region  = seq[q_start:q_end]

            # find BinI sites
            for pattern, site_name in [
                (BINI,    "BinI_CCTAGG"),
                (XAPI,    "XapI_AAATTCC"),
                (XAPI_RC, "XapI_GGAATTT"),
            ]:
                plen = len(pattern)
                for i in range(len(region) - plen + 1):
                    window = region[i:i+plen]
                    mismatches = sum(a != b for a, b in zip(window, pattern))
                    if mismatches <= 1:
                        # query position of this site
                        q_site_start = q_start + i
                        q_site_end   = q_start + i + plen - 1

                        # map back to reference
                        ref_site_start = query_to_ref.get(q_site_start)
                        ref_site_end   = query_to_ref.get(q_site_end)

                        if ref_site_start is None or ref_site_end is None:
                            continue

                        is_exact = mismatches == 0
                        site_rows.append({
                            "contig":      ru["contig"],
                            "ref_start":   min(ref_site_start, ref_site_end),
                            "ref_end":     max(ref_site_start, ref_site_end) + 1,
                            "site_name":   site_name,
                            "RU":          ru_idx,
                            "exact":       is_exact,
                            "mismatches":  mismatches,
                            "read_id":     read.query_name,
                        })

    bam.close()

    if not site_rows:
        return pd.DataFrame()

    df = pd.DataFrame(site_rows)

    # consensus: for each RU+site_name, take the most common reference position
    consensus = (
        df.groupby(["RU", "site_name", "ref_start", "ref_end", "contig", "exact"])
        .size()
        .reset_index(name="n_reads")
        .sort_values(["RU", "site_name", "n_reads"], ascending=[True, True, False])
        .drop_duplicates(subset=["RU", "site_name"], keep="first")
    )

    return consensus

def find_reference_in_allele_dir(allele_dir, allele_id):
    """Find reference FASTA in allele folder — Mode A custom ref or Mode B helper copy."""
    # Mode A — custom reference named after allele
    for ext in [".fa", ".fasta"]:
        ref = allele_dir / f"{allele_id}{ext}"
        if ref.exists():
            return ref

    # Mode B — helper reference copied into folder
    for ext in [".fa", ".fasta"]:
        candidates = [
            f for f in allele_dir.glob(f"*{ext}")
            if not f.name.startswith(allele_id)
            and ".raw." not in f.name
        ]
        if candidates:
            return candidates[0]

    return None
  
def run_single_allele(allele_dir, allele_id, out_dir):
    """Shared logic for processing one allele — used by both modes."""

    # find aligned BAM — exclude T2T subset and read_ids extraction BAMs
    exact = allele_dir / f"{allele_id}.bam"
    if exact.exists():
        bam_file = exact
    else:
        candidates = [
            f for f in allele_dir.glob("*.bam")
            if not f.name.endswith(".bam.bai")
            and "subset" not in f.name
            and "read_ids" not in f.name    # ← this was missing
        ]
        if not candidates:
            print(f"[WARN] No aligned BAM found in {allele_dir}, skipping.")
            return None, None
        # prefer BAM aligned to helper/custom reference over others
        ref_aligned = [f for f in candidates if f.stem.startswith(allele_id + "_")]
        bam_file = ref_aligned[0] if ref_aligned else candidates[0]

    # find curated TSV
    curated_tsv = allele_dir / f"{allele_id}.curated.tsv"
    if not curated_tsv.exists():
        candidates = list(allele_dir.glob("*.curated.tsv"))
        if not candidates:
            print(f"[WARN] No curated TSV found in {allele_dir}, skipping.")
            return None, None
        curated_tsv = candidates[0]

    print(f"[INFO] Processing: {bam_file.name} with {curated_tsv.name}")

    per_ru_df = process_allele(
        bam_path         = bam_file,
        curated_tsv_path = curated_tsv,
        out_dir          = out_dir,
        allele_id        = allele_id,
    )

    if per_ru_df is None or per_ru_df.empty:
        return None, None

    consensus_df = build_allele_bx_consensus(per_ru_df, allele_id)
    array_summary = None

    if not consensus_df.empty:
        consensus_df.to_csv(
            out_dir / f"{allele_id}_BX_consensus.csv",
            sep=";", index=False, encoding="utf-8-sig"
        )
        array_summary = build_allele_array_summary(consensus_df, allele_id)
        print(f"[INFO] Array structure: {array_summary.get('array_structure', '')}")

        # restriction site BED
        ru_intervals = get_ru_intervals_from_curated_tsv(curated_tsv)
        site_df = find_restriction_sites_in_ru(
            bam_path     = bam_file,
            ru_intervals = ru_intervals,
            prefix_name  = allele_id,
        )

        site_bed_path = out_dir / f"{allele_id}_BX_restriction_sites.bed"
        site_color_map = {
            "BinI_CCTAGG":  ("0,100,255", 4),
            "XapI_AAATTCC": ("255,50,0",  4),
            "XapI_GGAATTT": ("200,50,0",  4),
        }

        with open(site_bed_path, "w") as site_f:
            site_f.write(
                "track type=bed itemRgb=On name='BX_restriction_sites' "
                "description='BinI/XapI restriction site positions'\n"
            )
            if not site_df.empty:
                for _, srow in site_df.sort_values(["RU", "ref_start"]).iterrows():
                    contig    = str(srow["contig"])
                    s_start   = int(srow["ref_start"])
                    s_end     = int(srow["ref_end"])
                    site_name = str(srow["site_name"])
                    ru_num    = int(srow["RU"])
                    n_reads   = int(srow["n_reads"])
                    exact     = bool(srow["exact"])

                    label = (
                        f"RU{ru_num:02d}_{site_name}_"
                        f"{'exact' if exact else 'fuzzy'}(n={n_reads})"
                    )
                    color, score = site_color_map.get(site_name, ("128,128,128", 1))
                    if not exact:
                        color = "180,180,180"

                    site_f.write(
                        f"{contig}\t{s_start}\t{s_end}\t{label}\t"
                        f"{score}\t+\t{s_start}\t{s_end}\t{color}\n"
                    )

        print(f"[INFO] Restriction site BED: {site_bed_path}")

    return per_ru_df, array_summary

  
def main():
    ap = argparse.ArgumentParser(
        description="Per-RU BinI/XapI restriction site classification for DUCKS4 alleles."
    )

    # pipeline mode
    ap.add_argument("--bam-dir",    required=False, default="",
        help="Directory containing allele subfolders (allele_methylation/)")
    ap.add_argument("--manifest",   required=False, default="",
        help="TSV with columns allele_id and blast_file")

    # standalone mode
    ap.add_argument("--allele-dir", required=False, default="",
        help="Single finished ID2bam2meth output folder for standalone use")
    ap.add_argument("--allele-id",  required=False, default="",
        help="Allele ID for standalone mode (default: folder name)")

    ap.add_argument("--out-dir",    required=False, default="",
        help="Output directory for BX check results")

    args = ap.parse_args()

    # Standalone mode: single allele folder
    if args.allele_dir:
        allele_dir = Path(args.allele_dir)
        allele_id  = args.allele_id if args.allele_id else allele_dir.name
        out_dir    = Path(args.out_dir) if args.out_dir else allele_dir / "D4Z4_BX_check"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"[INFO] Standalone mode: {allele_id}")

        per_ru_df, array_summary = run_single_allele(allele_dir, allele_id, out_dir)

        if array_summary:
            pd.DataFrame([array_summary]).to_csv(
                out_dir / f"{allele_id}_array_structure.csv",
                sep=";", index=False, encoding="utf-8-sig"
            )

        print(f"[INFO] BX check complete: {out_dir}")
        return

    # Pipeline mode: multiple alleles via manifest
    if not args.bam_dir or not args.manifest:
        print("[ERROR] Provide either --allele-dir (standalone) or --bam-dir + --manifest (pipeline).")
        sys.exit(1)

    bam_dir = Path(args.bam_dir)
    out_dir = Path(args.out_dir) / "D4Z4_BX_check" if args.out_dir else bam_dir / "D4Z4_BX_check"
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest, sep="\t", dtype=str).fillna("")
    allele_blast_map = dict(zip(manifest["allele_id"], manifest["blast_file"]))

    all_summaries       = []
    all_array_summaries = []

    for allele in allele_blast_map.keys():
        allele_dir = bam_dir / allele

        if not allele_dir.is_dir():
            print(f"[WARN] No subfolder found for {allele}, skipping.")
            continue

        per_ru_df, array_summary = run_single_allele(allele_dir, allele, out_dir)

        if per_ru_df is not None:
            all_summaries.append(per_ru_df)
        if array_summary is not None:
            all_array_summaries.append(array_summary)

    if all_summaries:
        combined = pd.concat(all_summaries, ignore_index=True)
        combined.to_csv(
            out_dir / "all_alleles_BX_per_RU.csv",
            sep=";", index=False, encoding="utf-8-sig"
        )

    if all_array_summaries:
        array_df = pd.DataFrame(all_array_summaries)
        array_df.to_csv(
            out_dir / "all_alleles_array_structure.csv",
            sep=";", index=False, encoding="utf-8-sig"
        )
        print(f"[INFO] Array structure summary: {out_dir / 'all_alleles_array_structure.csv'}")

    if not all_summaries:
        print("[WARN] No results produced — check BAM and curated TSV paths.")

    print(f"[INFO] BX check complete: {out_dir}")


if __name__ == "__main__":
    main()
