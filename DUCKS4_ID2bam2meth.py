#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
import os
from pathlib import Path

import pandas as pd
import pysam


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def get_read_sequence(bam_path, read_id):
    bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
    seq = None
    is_reverse = False
    read_len = None
    is_unmapped = False

    for aln in bam.fetch(until_eof=True):
        if aln.query_name != read_id:
            continue

        if not aln.is_secondary and not aln.is_supplementary:
            seq = aln.query_sequence
            is_reverse = aln.is_reverse
            read_len = aln.query_length
            is_unmapped = aln.is_unmapped
            break

        if seq is None:
            seq = aln.query_sequence
            is_reverse = aln.is_reverse
            read_len = aln.query_length
            is_unmapped = aln.is_unmapped

    bam.close()

    if seq is None:
        raise ValueError(f"Read ID '{read_id}' not found in {bam_path}")

    return seq, is_reverse, read_len, is_unmapped


def write_fasta(read_id, sequence, fasta_path):
    with open(fasta_path, "w") as fh:
        fh.write(f">{read_id}\n")
        for i in range(0, len(sequence), 60):
            fh.write(sequence[i:i + 60] + "\n")


def index_fasta_with_samtools(fasta_path):
    try:
        subprocess.run(["samtools", "faidx", str(fasta_path)], check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: samtools faidx failed: {e}", file=sys.stderr)
        sys.exit(1)


def run_blast_on_fasta(fasta_path, blast_db, blast_out, blastn_bin, threads="32"):
    cmd = [
        str(blastn_bin),
        "-db", str(blast_db),
        "-query", str(fasta_path),
        "-num_threads", str(threads),
        "-out", str(blast_out),
        "-outfmt", "6 qseqid sseqid pident slen length mismatch gapopen qstart qend sstart send evalue bitscore",
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: BLAST failed: {e}", file=sys.stderr)
        sys.exit(1)

def extract_anchor_haplotype(name: str):
    s = str(name)
    for tag in ["4qA", "4qB", "10qA"]:
        if s.startswith(f"{tag}_D4F104S1") or s.startswith(f"{tag}_pLAM"):
            return tag
    return None


def classify_blast_subject(sseqid: str):
    s = str(sseqid)

    if "CLUHP4" in s:
        return "CLUHP4"
    if "DUX4_end" in s:
        return "DUX4_end"
    if "D4F104S1" in s:
        return "D4F104S1"
    if "pLAM" in s or "PLAM" in s:
        return "pLAM"

    # real repeat-like hits in fresh BLAST
    # examples: c4, c10, c4S, c10S, RU01_c4, 4q35_D4Z4, 10q26_D4Z4
    if s in {"c4", "c10", "c4S", "c10S"}:
        return "repeat"
    if s.startswith("RU"):
        return "repeat"
    if "D4Z4" in s:
        return "repeat"

    # pseudo/control-ish hits should not count as true repeats
    if "chr4_ctrl" in s or "chr10_ctrl" in s:
        return "ctrl"

    return "other"
  
def infer_orientation_from_blast(blast_path, target_read_id):
    """
    Decide orientation from the best framed block in fresh BLAST output.

    Preferred pattern:
      forward  = CLUHP4 ... D4F104S1 ... repeats ... pLAM ... DUX4_end
      inverted = DUX4_end ... pLAM ... repeats ... D4F104S1 ... CLUHP4
    """
    hits = []

    with open(blast_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            if len(cols) < 9:
                continue

            qseqid = cols[0].strip()
            if qseqid != target_read_id:
                continue

            sseqid = cols[1].strip()

            try:
                pident = float(cols[2])
                qstart = int(cols[7])
                qend = int(cols[8])
            except ValueError:
                continue

            start = min(qstart, qend)
            end = max(qstart, qend)

            hits.append({
                "name": sseqid,
                "class": classify_blast_subject(sseqid),
                "hap": extract_anchor_haplotype(sseqid),
                "pident": pident,
                "start": start,
                "end": end,
                "mid": (start + end) / 2.0,
            })

    if not hits:
        return "unknown"

    hit_df = pd.DataFrame(hits)

    marker_df = hit_df[hit_df["class"] == "D4F104S1"].copy()
    plam_df = hit_df[hit_df["class"] == "pLAM"].copy()
    repeat_df = hit_df[hit_df["class"] == "repeat"].copy()
    cluhp4_df = hit_df[hit_df["class"] == "CLUHP4"].copy()
    dux4end_df = hit_df[hit_df["class"] == "DUX4_end"].copy()

    if marker_df.empty or plam_df.empty:
        return "unknown"

    candidates = []

    for _, m in marker_df.iterrows():
        for _, p in plam_df.iterrows():

            if m["start"] < p["start"]:
                orientation = "forward"
                array_start = int(m["start"])
                array_end = int(p["end"])
                proximal_anchor_pos = int(m["start"])
                distal_anchor_pos = int(p["end"])
            else:
                orientation = "inverted"
                array_start = int(p["start"])
                array_end = int(m["end"])
                proximal_anchor_pos = int(m["end"])
                distal_anchor_pos = int(p["start"])

            span = array_end - array_start
            if span <= 0:
                continue

            # repeat hits inside the candidate array block
            repeats_between = repeat_df[
                (repeat_df["start"] >= array_start) &
                (repeat_df["end"] <= array_end)
            ].copy()

            n_repeats_between = len(repeats_between)

            # closest CLUHP4 and DUX4_end
            if not cluhp4_df.empty:
                dist_cluhp4_to_marker = min(abs(int(x) - proximal_anchor_pos) for x in cluhp4_df["mid"])
            else:
                dist_cluhp4_to_marker = 10**9

            if not dux4end_df.empty:
                dist_plam_to_dux4end = min(abs(int(x) - distal_anchor_pos) for x in dux4end_df["mid"])
            else:
                dist_plam_to_dux4end = 10**9

            # framing presence on correct sides
            if orientation == "forward":
                has_cluhp4_proximal = int(
                    not cluhp4_df.empty and (cluhp4_df["end"] <= m["start"]).any()
                )
                has_dux4end_distal = int(
                    not dux4end_df.empty and (dux4end_df["start"] >= p["end"]).any()
                )
            else:
                has_cluhp4_proximal = int(
                    not cluhp4_df.empty and (cluhp4_df["start"] >= m["end"]).any()
                )
                has_dux4end_distal = int(
                    not dux4end_df.empty and (dux4end_df["end"] <= p["start"]).any()
                )

            # soft haplotype bonus only
            hap_match = int(
                pd.notna(m["hap"]) and pd.notna(p["hap"]) and str(m["hap"]) == str(p["hap"])
            )

            mean_anchor_pid = (float(m["pident"]) + float(p["pident"])) / 2.0

            # penalties for huge distances from framing controls
            # pseudo hits far away (~40 kb) should lose strongly
            penalty_cluhp4 = dist_cluhp4_to_marker
            penalty_dux4end = dist_plam_to_dux4end

            # prefer plausible compact blocks, but do not hard-cut
            plausible_span = int(1000 <= span <= 250000)

            # main score: framing + repeats + distance structure
            score = 0
            score += has_cluhp4_proximal * 1000
            score += has_dux4end_distal * 1000
            score += min(n_repeats_between, 100) * 50
            score += plausible_span * 200
            score += hap_match * 25
            score += mean_anchor_pid

            score -= penalty_cluhp4 / 100.0
            score -= penalty_dux4end / 100.0

            # mild penalty for zero repeats between anchors
            if n_repeats_between == 0:
                score -= 500

            candidates.append({
                "orientation": orientation,
                "score": score,
                "span": span,
                "marker_name": m["name"],
                "plam_name": p["name"],
                "marker_start": int(m["start"]),
                "plam_start": int(p["start"]),
                "n_repeats_between": n_repeats_between,
                "dist_cluhp4_to_marker": dist_cluhp4_to_marker,
                "dist_plam_to_dux4end": dist_plam_to_dux4end,
                "has_cluhp4_proximal": has_cluhp4_proximal,
                "has_dux4end_distal": has_dux4end_distal,
                "hap_match": hap_match,
                "mean_anchor_pid": mean_anchor_pid,
            })

    if not candidates:
        return "unknown"

    cand_df = pd.DataFrame(candidates).sort_values(
        by=[
            "score",
            "has_cluhp4_proximal",
            "has_dux4end_distal",
            "n_repeats_between",
            "mean_anchor_pid",
        ],
        ascending=[False, False, False, False, False],
    )

    best = cand_df.iloc[0]

    print(
        "Best framed block for orientation: "
        f"{best['marker_name']} <-> {best['plam_name']} | "
        f"orientation={best['orientation']} | "
        f"repeats_between={best['n_repeats_between']} | "
        f"dist(CLUHP4,marker)={best['dist_cluhp4_to_marker']} | "
        f"dist(pLAM,DUX4_end)={best['dist_plam_to_dux4end']} | "
        f"score={best['score']:.2f}"
    )

    return best["orientation"]

def parse_blast_to_bed(blast_path, target_read_id, bed_path):
    """
    Parse final-orientation BLAST output into raw BED.
    """
    hits = {}

    with open(blast_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            cols = line.split()
            if len(cols) < 9:
                continue

            qseqid = cols[0].strip()
            if qseqid != target_read_id:
                continue

            region = cols[1].strip().rstrip(".")
            
            # label pseudo/low-identity pLAM separately
            if ("pLAM" in region or "PLAM" in region) and pid < 90:
                region = region.replace("_pLAM", "_pLAM_lowpid")
                region = region.replace("_PLAM", "_pLAM_lowpid")

            try:
                pid = float(cols[2].replace(",", "."))
                qstart = int(cols[7])
                qend = int(cols[8])
            except ValueError:
                continue

            start = min(qstart, qend) - 1
            end = max(qstart, qend)

            key = (start, end)
            best = hits.get(key)
            if best is None or pid > best["pid"]:
                hits[key] = {"pid": pid, "region": region}

    if not hits:
        raise ValueError(f"No BLAST hits found for read {target_read_id} in {blast_path}")

    with open(bed_path, "w") as out:
        for (start, end) in sorted(hits.keys()):
            out.write(f"{target_read_id}\t{start}\t{end}\t{hits[(start, end)]['region']}\n")


def run_annotation_curation(curate_script, raw_bed, curated_tsv, curated_bed):
    cmd = [
        sys.executable,
        str(curate_script),
        "--bed", str(raw_bed),
        "--out-curation", str(curated_tsv),
        "--out-bed", str(curated_bed),
        "--verbose",
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: annotation BED curation failed: {e}", file=sys.stderr)
        sys.exit(1)


def checkfile_fastq(file_in):
    file_out = file_in.rsplit('.', 1)[0] + '.fastq'

    if file_in.endswith('bam'):
        print(".bam-file is being converted to .fastq-file with samtools.")
        print("       ")
        try:
            with open(file_out, 'w') as out_f:
                result = subprocess.run(
                    ['samtools', 'fastq', '-T', '*', file_in],
                    check=True, text=True, stdout=out_f, stderr=subprocess.PIPE
                )
            print(f"The file was successfully converted to {file_out}.")
            if result.stderr:
                print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Error while converting the file: {e}")
            print(f"{e.stderr}")
            return None

    elif file_in.endswith('gz'):
        print(".fastq.gz-file is being unzipped.")
        print("       ")
        try:
            with open(file_out, 'w') as out_f:
                result = subprocess.run(
                    ['gunzip', '-c', file_in],
                    check=True, text=True, stdout=out_f, stderr=subprocess.PIPE
                )
            print(f"The file was successfully converted to {file_out}.")
            if result.stderr:
                print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
            print(f"Fehler bei der Umwandlung der Datei: {e}")
            print(f"Standardfehlerausgabe: {e.stderr}")
            return None

    elif file_in.endswith('fastq'):
        file_out = file_in

    else:
        print("This file-format is not supported!")
        return None

    return file_out


# combined region methylation helpers


def read_curated_annotation_tsv(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    for col in ["start", "end", "length_bp", "ru_index"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def read_modkit_stats_tsv(path):
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")

    rename_map = {
        "#chrom": "contig",
        "chrom": "contig",
        "name": "display_name",
        "count_m": "total_nmod",
        "count_valid_m": "total_valid_cov",
        "percent_m": "weighted_mean_methylation",
    }

    df = df.rename(columns=rename_map)

    if "contig" not in df.columns:
        first_col = df.columns[0]
        df = df.rename(columns={first_col: "contig"})

    required = ["contig", "start", "end", "weighted_mean_methylation"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(
            "modkit stats output is missing required columns after normalization: "
            + ", ".join(missing)
            + "\nFound columns: "
            + ", ".join(map(str, df.columns))
        )

    numeric_cols = [
        "start", "end",
        "total_nmod", "total_valid_cov", "weighted_mean_methylation",
        "count_h", "count_valid_h", "percent_h",
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    if "total_ncanonical" not in df.columns:
        if {"total_valid_cov", "total_nmod"}.issubset(df.columns):
            df["total_ncanonical"] = df["total_valid_cov"] - df["total_nmod"]
        else:
            df["total_ncanonical"] = pd.NA

    if "mean_methylation" not in df.columns and "weighted_mean_methylation" in df.columns:
        df["mean_methylation"] = df["weighted_mean_methylation"]

    if "n_sites" not in df.columns:
        df["n_sites"] = pd.NA

    if "median_methylation" not in df.columns:
        df["median_methylation"] = df["weighted_mean_methylation"]

    if "min_methylation" not in df.columns:
        df["min_methylation"] = df["weighted_mean_methylation"]

    if "max_methylation" not in df.columns:
        df["max_methylation"] = df["weighted_mean_methylation"]

    return df


def normalize_bool_series(series):
    return (
        series.astype(str)
        .str.strip()
        .str.upper()
        .isin(["TRUE", "WAHR", "1"])
    )


def merge_curated_with_modkit(curated_df, modkit_df):
    key_cols = ["contig", "start", "end"]
    merged = curated_df.merge(
        modkit_df,
        on=key_cols,
        how="left",
        suffixes=("", "_meth"),
    )
    return merged


# distal_unit is added to the curated regions BED before modkit stats,
# distal_unit: last-repeat + pLAM region with PAS inside for combined methylation call.

def main():
    script_path = os.path.dirname(os.path.realpath(__file__))
    db_path = os.path.join(script_path, "ressources", "blast_db", "FSHD-blast")
    blast_path = os.path.join(script_path, "ressources", "tools", "ncbi-blast-2.14.0+", "bin", "blastn")
    curate_script = os.path.join(script_path, "ressources", "tools", "curate_annotation_bed.py")

    ap = argparse.ArgumentParser(
        description="""
        DUCKS4 - ID2bam2meth

        Reference:
        A) create custom reference from a read:
           --id_ref read-id --bam_ref read.bam
        B) provide existing reference:
           --ref ref.fasta

        Alignment:
          --bam reads.bam
          [--txt] read_ids.txt

        Methylation:
          --methyl
            - Mode A: read-id: uses blast-annotation/curation for methylation
            - Mode B: ref.fasta: uses blast-annotation/curation for methylation if no --region/--regions_bed is given
            
        BinI/XapI - check (EXPERIMENTAL)
        New BinI/XapI check included, 
        classifies each D4Z4 RU as Chr4_D4Z4 (B-/X+), Chr10_D4Z4 (B+/X-), or Hybrid_D4Z4 (B-/X-)
        based on the presence of BinI/AvrII (CCTAGG) and XapI/ApoI (AAATTCC) restriction sites 
        detected directly from read alignments.
        Per-RU consensus is built by majority vote across all covering reads, producing BED files
        for IGV visualization (identified RUs and used Restriction-sites for evaluation) and summary CSVs.
        Attention: it is an experimental feature, please use and evaluate with care.
        --skip-bx     # switch off this feature
        """
    )

    ap.add_argument("--id_ref", required=False, help="Read-id used as custom reference (Mode A).")
    ap.add_argument("--bam_ref", required=False, help="BAM containing the reference-read (Mode A).")
    ap.add_argument("--ref", required=False, help="Use an existing reference FASTA (Mode B).")
    ap.add_argument("--bam", required=False, help="Input BAM to extract/align reads against the reference.")
    ap.add_argument("--txt", required=False, help="Optional read_id.txt to subset reads from --bam.")
    ap.add_argument("--methyl", required=False, action='store_true', help="Run modkit pileup + stats.")
    ap.add_argument("--region", required=False, help="Optional region. Format: chr:start-end")
    ap.add_argument("--regions_bed", required=False, help="Optional BED with multiple regions.")
    ap.add_argument("--skip-bx", required=False, action="store_true", help="Skip BX restriction site check.")
    ap.add_argument("--threads", dest="threads", default="45", required=False, help="Threads (default 45).")
    ap.add_argument("--out_prefix", required=False, default=None, help="Prefix name (default: id_ref or ref basename).")
    ap.add_argument("--out_path", required=False, help="Output directory. Default: created next to BAM/ref.")

    args = ap.parse_args()

    print("       ")
    print("####   DUCKS4 - ID2bam2meth   ####")
    print("Version 1.1.0")
    print("       ")

    modeA = bool(args.id_ref and args.bam_ref)
    modeB = bool(args.ref)

    if not modeA and not modeB:
        print("ERROR: Choose either:")
        print("  Mode A: --id_ref --bam_ref")
        print("  Mode B: --ref")
        sys.exit(1)

    if modeA and modeB:
        print("ERROR: Please choose either Mode A or Mode B, not both.")
        sys.exit(1)

    if args.txt and not args.bam:
        print("ERROR: --bam is required when --txt is provided.")
        sys.exit(1)

    if args.bam is None and args.methyl:
        print("ERROR: --bam is required if you want to run --methyl.")
        sys.exit(1)

    if args.out_prefix:
        prefix_name = args.out_prefix
    else:
        prefix_name = args.id_ref if modeA else Path(args.ref).stem

    if args.out_path:
        outpath = Path(args.out_path)
        outpath.mkdir(parents=True, exist_ok=True)
    else:
        if args.bam:
            base_dir = Path(args.bam).resolve().parent
        elif modeA:
            base_dir = Path(args.bam_ref).resolve().parent
        else:
            base_dir = Path(args.ref).resolve().parent

        dir_name = f"custref2bam_{prefix_name}" if modeA else f"ref2bam_{prefix_name}"
        outpath = base_dir / dir_name
        outpath.mkdir(parents=True, exist_ok=True)

    def minimap2(file, ref_fasta, out_dir):
        print("       ")
        print("Alignment with Minimap2:")
        print("       ")
        file_in = os.path.basename(file)
        ref_name = os.path.basename(ref_fasta).split('.')[0]
        if args.ref:
            sam_file = f"{prefix_name}_{ref_name}.sam"
        else:
            sam_file = f"{prefix_name}.sam"

        subprocess.call([
            "minimap2", "-ax", "lr:hq", "--MD", "-L",
            "-t", str(args.threads), "-Y", "-y",
            "-o", str(Path(out_dir) / sam_file),
            str(ref_fasta), str(Path(file).resolve())
        ])
        return sam_file

    def samtools_bam(sam_input):
        print("       ")
        print("Sorting and indexing of .sam-file:")
        print("       ")
        sam_path = Path(sam_input)
        bam_file = sam_path.with_suffix(".bam").name
        bam_path = sam_path.with_suffix(".bam")

        subprocess.call(["samtools", "sort", "-m", "15G", "-o", str(bam_path), str(sam_path)])
        subprocess.call(["samtools", "index", str(bam_path)])
        subprocess.call(["rm", str(sam_path)])
        return bam_file

    ref_fasta = None
    bed_path = None
    final_regions_bed = None
    curated_tsv_path = None
    curated_bed_path = None

    # Mode A
    if modeA:
        bam_ref = Path(args.bam_ref).resolve()

        raw_fasta_path = outpath / f"{prefix_name}.raw.fa"
        final_fasta_path = outpath / f"{prefix_name}.fa"
        bed_path = outpath / f"{prefix_name}.bed"

        curated_tsv_path = outpath / f"{prefix_name}.curated.tsv"
        curated_bed_path = outpath / f"{prefix_name}.curated.bed"

        raw_blast_path = outpath / f"{prefix_name}.raw.blast.tsv"
        final_blast_path = outpath / f"{prefix_name}.final.blast.tsv"

        print("Create custom reference from read:")
        print(f"  id_ref    : {args.id_ref}")
        print(f"  bam_ref   : {bam_ref}")
        print(f"  blast_db  : {db_path}")
        print(f"  blast_bin : {blast_path}")
        print(f"  curate    : {curate_script}")

        seq, is_reverse, read_len, is_unmapped = get_read_sequence(bam_ref, args.id_ref)

        print(f"BAM flag info: is_reverse={is_reverse}, is_unmapped={is_unmapped}")
        print("Reference orientation will be decided from re-BLAST annotation order.")

        write_fasta(args.id_ref, seq, raw_fasta_path)
        index_fasta_with_samtools(raw_fasta_path)

        run_blast_on_fasta(
            fasta_path=raw_fasta_path,
            blast_db=db_path,
            blast_out=raw_blast_path,
            blastn_bin=blast_path,
            threads=args.threads,
        )

        orientation = infer_orientation_from_blast(raw_blast_path, args.id_ref)
        print(f"Orientation inferred from framed BLAST block: {orientation}")

        flipped = False
        if orientation == "inverted":
            print("Re-BLAST anchor order indicates inverted orientation: pLAM ... D4F104S1")
            print("Reverse-complementing read to standardize reference to D4F104S1 -> D4Z4 -> pLAM.")
            seq = reverse_complement(seq)
            flipped = True
        elif orientation == "forward":
            print("Re-BLAST anchor order indicates forward orientation: D4F104S1 ... pLAM")
            print("Keeping read sequence as it is.")
        else:
            print("WARNING: Could not infer anchor orientation from re-BLAST.")
            print("Keeping read sequence as it is.")

        write_fasta(args.id_ref, seq, final_fasta_path)
        index_fasta_with_samtools(final_fasta_path)

        run_blast_on_fasta(
            fasta_path=final_fasta_path,
            blast_db=db_path,
            blast_out=final_blast_path,
            blastn_bin=blast_path,
            threads=args.threads,
        )

        parse_blast_to_bed(final_blast_path, args.id_ref, bed_path)

        if not Path(curate_script).exists():
            print(f"ERROR: curate script not found: {curate_script}", file=sys.stderr)
            sys.exit(1)

        print("Running annotation BED curation.")
        run_annotation_curation(
            curate_script=curate_script,
            raw_bed=bed_path,
            curated_tsv=curated_tsv_path,
            curated_bed=curated_bed_path,
        )

        ref_fasta = final_fasta_path
        final_regions_bed = curated_bed_path

        for fp in [
            raw_fasta_path,
            Path(str(raw_fasta_path) + ".fai"),
            raw_blast_path,
            bed_path,
        ]:
            try:
                if Path(fp).exists():
                    Path(fp).unlink()
            except Exception as e:
                print(f"WARNING: could not delete temporary file {fp}: {e}")

        print("Done creating custom reference.")
        print(f"FINAL FASTA       : {final_fasta_path}")
        print(f"FINAL FAIDX       : {final_fasta_path}.fai")
        print(f"FINAL BLAST       : {final_blast_path}")
        print(f"CURATED TSV       : {curated_tsv_path}")
        print(f"CURATED BED       : {curated_bed_path}")
        print(f"FINAL ORIENT      : {'flipped_to_forward' if flipped else 'kept_as_is'}")

    # Mode B
    else:
        ref_fasta = Path(args.ref).resolve()
        ref_name = ref_fasta.stem
        print("Using provided reference FASTA:")
        print(f"  REF : {ref_fasta}")
        
        # copy reference into output folder for traceability 
        import shutil
        ref_copy = outpath / ref_fasta.name
        if not ref_copy.exists() or ref_copy.stat().st_size < 1000:
            try:
                shutil.copy2(str(ref_fasta), str(ref_copy))
                print(f"  Reference copied to output: {ref_copy.name}")
            except Exception as e:
                print(f"  [WARN] Could not copy reference: {e}")
        ref_fai = Path(str(ref_fasta) + ".fai")
        if ref_fai.exists():
            fai_copy = outpath / ref_fai.name
            if not fai_copy.exists():
                try:
                    shutil.copy2(str(ref_fai), str(fai_copy))
                except Exception as e:
                    print(f"  [WARN] Could not copy reference index: {e}")
        
        def get_first_fasta_header(fasta_path):
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        return line[1:].strip().split()[0]
            raise ValueError(f"No FASTA header found in {fasta_path}")
          
        ref_read_id = get_first_fasta_header(ref_fasta)
    
        if not Path(str(ref_fasta) + ".fai").exists():
            print("No .fai found -> indexing reference with samtools faidx.")
            index_fasta_with_samtools(ref_fasta)
    
        # optional BLAST-based annotation for provided reference FASTA
        try:
    
            bed_path = outpath / f"{prefix_name}_{ref_name}.bed"
            curated_tsv_path = outpath / f"{prefix_name}_{ref_name}.curated.tsv"
            curated_bed_path = outpath / f"{prefix_name}_{ref_name}.curated.bed"
            final_blast_path = outpath / f"{prefix_name}_{ref_name}.final.blast.tsv"
    
            print("BLASTing provided reference FASTA for automatic annotation.")
            run_blast_on_fasta(
                fasta_path=ref_fasta,
                blast_db=db_path,
                blast_out=final_blast_path,
                blastn_bin=blast_path,
                threads=args.threads,
            )
    
            parse_blast_to_bed(final_blast_path, ref_read_id, bed_path)
    
            run_annotation_curation(
                curate_script=curate_script,
                raw_bed=bed_path,
                curated_tsv=curated_tsv_path,
                curated_bed=curated_bed_path,
            )
    
            final_regions_bed = curated_bed_path
    
            print(f"CURATED TSV       : {curated_tsv_path}")
            print(f"CURATED BED       : {curated_bed_path}")
    
        except Exception as e:
            print(f"WARNING: automatic BLAST/curation for provided reference failed: {e}")
            final_regions_bed = None
            
    if args.methyl and modeB and (not args.regions_bed) and final_regions_bed is None and (not args.region):
        print("ERROR: Could not determine regions for --ref + --methyl.")
        print("Please provide --region or --regions_bed, or ensure automatic BLAST annotation works.")
        sys.exit(1)

    if args.bam:
        bam_in = Path(args.bam).resolve()

        if args.txt:
            print("Extract reads from bam using read_id.txt.")
            txt_name = Path(args.txt).stem
            bam_name = bam_in.stem
            ID_bam = f"{bam_name}_{txt_name}.bam"
            ID_out = outpath / ID_bam

            try:
                with open(ID_out, 'w') as out_f:
                    subprocess.run(
                        ['samtools', 'view', '-h', '-N', args.txt, str(bam_in)],
                        check=True, text=True, stdout=out_f, stderr=subprocess.PIPE
                    )
                print(f"Extracted BAM written to: {ID_out}")
            except subprocess.CalledProcessError as e:
                print(f"Error with samtools view: {e}")
                sys.exit(1)

            map_input = str(ID_out)
        else:
            print("No --txt provided -> whole bam will be aligned against reference.")
            map_input = str(bam_in)

        print("Align reads to reference.")
        fastq = checkfile_fastq(map_input)
        if fastq is None:
            print("ERROR: could not create fastq.")
            sys.exit(1)

        sam_map = minimap2(fastq, ref_fasta, outpath)
        bam_map = samtools_bam(str(outpath / sam_map))

        if map_input.endswith(".bam") and fastq.endswith(".fastq") and Path(fastq).exists():
            subprocess.call(["rm", fastq])

        print(f"Aligned BAM: {outpath / bam_map}")

    else:
        bam_map = None
        print("No alignment requested (no --bam).")

    if args.methyl:
        if bam_map is None:
            print("ERROR: --methyl requested but no aligned BAM exists.")
            sys.exit(1)

        print("       ")
        print("Start Methylation-analysis with Modkit.")
        print("       ")

        bam_map_name = Path(bam_map).stem
        meth_path = outpath / f"{prefix_name}_methylation-analysis"
        meth_path.mkdir(parents=True, exist_ok=True)
        
        if args.ref:
            modkit_bed = f"{prefix_name}_{ref_name}_modkit-methyl.bed"
        else:
            modkit_bed = f"{prefix_name}_modkit-methyl.bed"
        modkit_bed_path = meth_path / modkit_bed

        subprocess.call([
            "modkit", "pileup",
            str(outpath / bam_map),
            str(modkit_bed_path),
            "--cpg",
            "--ref", str(ref_fasta)
        ])

        subprocess.call(["bgzip", str(modkit_bed_path)])
        subprocess.call(["tabix", str(modkit_bed_path) + ".gz"])
        
        if args.regions_bed:
            regions_bed = Path(args.regions_bed).resolve()
        elif final_regions_bed is not None and Path(final_regions_bed).exists():
            regions_bed = final_regions_bed
        else:
            chrom, coords = args.region.split(":")
            start, end = map(int, coords.split("-"))
            regions_bed = meth_path / "coordinates_methcalc.bed"
            with open(regions_bed, "w") as f:
                f.write(f"{chrom}\t{start - 1}\t{end}\n")

        modout_all = meth_path / "modkit-STATS.tsv"
        subprocess.call([
            "modkit", "stats",
            "--regions", str(regions_bed),
            "-o", str(modout_all),
            str(modkit_bed_path) + ".gz"
        ])



        modkit_df = read_modkit_stats_tsv(modout_all)
        
        if args.ref:
            output_bg = f"{prefix_name}_{ref_name}.bedgraph"
            output_bed = f"{prefix_name}_{ref_name}.bed"
        else:
            output_bg = f"{prefix_name}.bedgraph"
            output_bed = f"{prefix_name}.bed"

        with open(meth_path / output_bg, "w") as fout:
            fout.write("track type=bedGraph name='Percent_Methylation' description='Percent Methylation'\n")
            for _, row in modkit_df.iterrows():
                if pd.isna(row.get("weighted_mean_methylation")):
                    continue
                fout.write(
                    f"{row['contig']}\t{int(row['start'])}\t{int(row['end'])}\t{float(row['weighted_mean_methylation']):.6f}\n"
                )

        with open(meth_path / output_bed, "w") as fout:
            fout.write("track type=bed name='Methylation_Labels' description='Percent Methylation as Name'\n")
            for _, row in modkit_df.iterrows():
                if pd.isna(row.get("weighted_mean_methylation")):
                    continue
                fout.write(
                    f"{row['contig']}\t{int(row['start'])}\t{int(row['end'])}\t{float(row['weighted_mean_methylation']):.2f}\n"
                )


        subprocess.call(["chown", "-R", "777", str(meth_path)])

        print("Methylation outputs:")
        print(f"  {modout_all}")
        print(f"  {meth_path / output_bg}")
        print(f"  {meth_path / output_bed}")

        
    # BinI/XapI per-RU check (experimental)
    bx_script = os.path.join(script_path, "ressources", "tools", "BX_check.py")

    if (
          not args.skip_bx
          and os.path.exists(bx_script)
          and curated_tsv_path is not None
          and Path(curated_tsv_path).exists()
        ):
          print("       ")
          print("Running BinI/XapI per-RU restriction site check.")

          bx_out_dir = outpath / "D4Z4_BX_check"
          bx_out_dir.mkdir(parents=True, exist_ok=True)

          try:
              subprocess.run([
                  "python3", bx_script,
                  "--allele-dir", str(outpath),
                  "--allele-id",  prefix_name,
                  "--out-dir",    str(bx_out_dir),
              ], check=True)

              # write BX site BED for IGV
              consensus_csv = bx_out_dir / f"{prefix_name}_BX_consensus.csv"

              if consensus_csv.exists():
                  curated_df = pd.read_csv(curated_tsv_path, sep="\t", dtype=str).fillna("")
                  for col in ["start", "end", "ru_index"]:
                      curated_df[col] = pd.to_numeric(curated_df[col], errors="coerce")

                  ru_df = curated_df[
                      curated_df["count_for_RU"].astype(str).str.upper().isin(["TRUE", "1"])
                  ].copy()

                  consensus_df = pd.read_csv(consensus_csv, sep=";", dtype=str).fillna("")
                  consensus_df["RU"] = pd.to_numeric(consensus_df["RU"], errors="coerce")

                  ru_df["ru_index"] = pd.to_numeric(ru_df["ru_index"], errors="coerce")
                  merged = consensus_df.merge(
                      ru_df[["ru_index", "start", "end", "contig"]],
                      left_on="RU", right_on="ru_index",
                      how="left"
                  )

                  bx_bed_path = bx_out_dir / f"{prefix_name}_BX_sites.bed"
                  bx_bg_path  = bx_out_dir / f"{prefix_name}_BX_sites.bedgraph"

                  color_map = {
                      "Chr10_D4Z4 (B+/X-)": "0,0,255",
                      "Chr4_D4Z4 (B-/X+)":  "255,0,0",
                      "Hybrid_D4Z4 (B-/X-)": "128,0,128",
                      "ambiguous (B+/X+)":  "128,128,128",
                      "unclassified":       "200,200,200",
                      "D4Z4-S":             "0,160,80",
                  }
                  
                  type_score = {
                      "Chr10_D4Z4 (B+/X-)": 3,
                      "Chr4_D4Z4 (B-/X+)":  2,
                      "Hybrid_D4Z4 (B-/X-)": 1,
                      "ambiguous (B+/X+)":  0,
                      "unclassified":       0,
                      "D4Z4-S":             2,
                  }

                  with open(bx_bed_path, "w") as bed_f, open(bx_bg_path, "w") as bg_f:
                      bed_f.write("track type=bed itemRgb=On name='BX_D4Z4_type' description='BinI/XapI D4Z4 unit classification'\n")
                      bg_f.write("track type=bedGraph name='BX_confidence' description='BX consensus confidence'\n")

                      for _, row in merged.sort_values("RU").iterrows():
                          if pd.isna(row.get("start")) or pd.isna(row.get("end")):
                              continue

                          contig     = str(row.get("contig", prefix_name))
                          start      = int(row["start"])
                          end        = int(row["end"])
                          d4z4_type  = str(row.get("consensus_type", "unclassified"))
                          b_status   = str(row.get("B_consensus", "?"))
                          x_status   = str(row.get("X_consensus", "?"))
                          confidence = row.get("confidence", "")
                          ru_num     = int(row["RU"]) if pd.notna(row["RU"]) else 0
                          flag       = str(row.get("flag", ""))

                          type_short = {
                              "Chr10_D4Z4 (B+/X-)": "c10",
                              "Chr4_D4Z4 (B-/X+)":  "c4",
                              "Hybrid_D4Z4 (B-/X-)": "hybrid",
                              "ambiguous (B+/X+)":  "ambiguous",
                              "unclassified":       "?",
                              "D4Z4-S":             "D4Z4-S",
                          }
                          d4z4_short = type_short.get(d4z4_type, "?")
                          label = f"RU{ru_num:02d}_{b_status}/{x_status}_{d4z4_short}"
                          if flag:
                              label += f"_({flag})"

                          color = color_map.get(d4z4_type, "200,200,200")
                          score = type_score.get(d4z4_type, 0)

                          bed_f.write(
                              f"{contig}\t{start}\t{end}\t{label}\t"
                              f"{score}\t+\t{start}\t{end}\t{color}\n"
                          )

                          if confidence != "" and not pd.isna(confidence):
                              try:
                                  bg_f.write(f"{contig}\t{start}\t{end}\t{float(confidence):.3f}\n")
                              except (ValueError, TypeError):
                                  pass

                  print(f"  BX site BED:            {bx_bed_path}")
                  print(f"  BX confidence bedGraph: {bx_bg_path}")

          except Exception as e:
              print(f"WARNING: BX check failed: {e}")          
    
    subprocess.call(["chmod", "-R", "777", outpath])

    print("       ")
    print("Workflow finished. Thank you for using this pipeline. TL :-)")
    print("       ")


if __name__ == "__main__":
    main()
