#!/usr/bin/env python3
import argparse, os, re
import pysam
import pandas as pd

def is_snv(rec) -> bool:
    """True if all ALT alleles are single-base and REF is single-base."""
    if rec.alts is None:
        return False
    if len(rec.ref) != 1:
        return False
    return all(len(alt) == 1 for alt in rec.alts)

def looks_fshd_annot(rec, ann_key="ANN", info_key="FSHD", ann_regex=r"FSHD"):
    """
    Heuristics:
      - INFO contains 'FSHD' (flag or value), OR
      - any entry in ANN contains 'FSHD' (case-insensitive).
    Adjust keys/regex as needed for your pipeline.
    """
    # 1) dedicated INFO tag (flag or value)
    if info_key in rec.info:
        val = rec.info[info_key]
        try:
            # INFO flag becomes True/False; others may be tuple/list/str
            if (isinstance(val, (bool, int)) and bool(val)) or (isinstance(val, (tuple, list)) and len(val) > 0) or (isinstance(val, str) and val.strip()):
                return True
        except Exception:
            return True  # present is enough

    # 2) search in ANN strings
    if ann_key in rec.info:
        rx = re.compile(ann_regex, re.IGNORECASE)
        ann = rec.info[ann_key]
        # ANN may be a tuple of strings; each is pipe-delimited by SnpEff
        for entry in ann:
            if rx.search(entry or ""):
                return True

    return False

def parse_ANN_best(ann_str):
    """
    Parse a single SnpEff ANN string (pipe-delimited). Returns a few fields.
    Format: Allele|Annotation|Impact|Gene_Name|Gene_ID|Feature_Type|Feature_ID|Transcript_BioType|Rank/Total|...
    """
    fields = (ann_str or "").split("|")
    out = {
        "ann_annotation": fields[1] if len(fields) > 1 else "",
        "ann_impact":     fields[2] if len(fields) > 2 else "",
        "ann_gene":       fields[3] if len(fields) > 3 else "",
        "ann_gene_id":    fields[4] if len(fields) > 4 else "",
        "ann_feature_id": fields[6] if len(fields) > 6 else "",
    }
    return out

def main():
    ap = argparse.ArgumentParser(description="Filter SNVs with FSHD annotation from a SnpEff/SnpSift VCF.")
    ap.add_argument("--vcf", required=True, help="Input VCF (can be .vcf or .vcf.gz)")
    ap.add_argument("--out-vcf", required=True, help="Output filtered VCF (.vcf or .vcf.gz)")
    ap.add_argument("--out-tsv", required=True, help="Output TSV table with annotations")
    ap.add_argument("--ann-key", default="ANN", help="INFO key for SnpEff annotations (default: ANN)")
    ap.add_argument("--info-key", default="FSHD", help="INFO key/flag marking FSHD hits (default: FSHD)")
    ap.add_argument("--ann-regex", default=r"FSHD", help="Regex to detect FSHD in ANN strings (default: FSHD)")
    args = ap.parse_args()

    infile = pysam.VariantFile(args.vcf)
    outfile = pysam.VariantFile(args.out_vcf, mode="w", header=infile.header)

    rows = []
    for rec in infile:
        if not is_snv(rec):
            continue
        if not looks_fshd_annot(rec, ann_key=args.ann_key, info_key=args.info_key, ann_regex=args.ann_regex):
            continue

        # write VCF record unchanged into filtered VCF
        outfile.write(rec)

        # build one or more table rows (one per ALT)
        ann_entries = []
        if args.ann_key in rec.info:
            ann_entries = list(rec.info[args.ann_key])  # tuple of strings
        # pick the first ANN as “best” for table (or aggregate if you want)
        best_ann = parse_ANN_best(ann_entries[0]) if ann_entries else {}

        base = {
            "chrom": rec.contig,
            "pos": rec.pos,                # 1-based
            "id": rec.id or "",
            "ref": rec.ref,
            "alt_all": ",".join(rec.alts or []),
            "qual": rec.qual if rec.qual is not None else "",
            "filter": ";".join(list(rec.filter.keys())) if rec.filter is not None else "PASS",
        }
        # add genotypes for all samples (GT only; extend if needed)
        for sample, call in rec.samples.items():
            gt = call.get("GT")
            base[f"{sample}_GT"] = "/".join("." if a is None else str(a) for a in (gt or []))

        # merge with parsed ANN fields
        base.update(best_ann)
        rows.append(base)

    infile.close()
    outfile.close()

    # Write table
    df = pd.DataFrame(rows)
    df.to_csv(args.out_tsv, sep="\t", index=False)

if __name__ == "__main__":
    main()
