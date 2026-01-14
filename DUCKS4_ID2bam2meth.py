#!/usr/bin/env python3
import argparse
import csv
import subprocess
import sys
import os, time, shutil
from pathlib import Path
from pprint import pprint

import pysam


def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def get_read_sequence(bam_path, read_id):
    bam = pysam.AlignmentFile(str(bam_path), "rb", check_sq=False)
    seq = None
    is_reverse = False
    read_len = None
    for aln in bam.fetch(until_eof=True):
        if aln.query_name != read_id:
            continue
        # Prefer primary, non-supplementary alignment
        if not aln.is_secondary and not aln.is_supplementary:
            seq = aln.query_sequence
            is_reverse = aln.is_reverse
            read_len = aln.query_length
            is_unmapped = aln.is_unmapped
            break
        # take first hit if no primary is found
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


def parse_blast_to_bed(blast_path, target_read_id, bed_path, read_len, flipped):
    hits = {}
    with open(blast_path) as fh:
        first = fh.readline().rstrip("\n")
        if not first:
            raise ValueError(f"BLAST file {blast_path} is empty")
        has_header = first.lower().startswith("read.id") or first.lower().startswith('"read.id"')
        # delimiter detection
        if "\t" in first:
            delim = "\t"
        elif ";" in first:
            delim = ";"
        elif "," in first:
            delim = ","
        else:
            delim = "\t"
        fh.seek(0)
        if has_header:
            # headered CSV/TSV 
            reader = csv.reader(fh, delimiter=delim)
            header = next(reader)
            header = [h.strip().strip('"') for h in header]
            required = ["read.id", "region", "percent.identity", "qstart", "qend"]
            missing = [c for c in required if c not in header]
            if missing:
                raise ValueError(
                    f"BLAST file missing required columns: {', '.join(missing)}\n"
                    f"Found columns: {header}"
                )
            idx = {name: header.index(name) for name in header}
            for cols in reader:
                if not cols:
                    continue
                # pad short lines
                if len(cols) < len(header):
                    cols = cols + [""] * (len(header) - len(cols))
                read_id = cols[idx["read.id"]].strip()
                if read_id != target_read_id:
                    continue
                region = cols[idx["region"]].strip().rstrip(".")
                pid_str = cols[idx["percent.identity"]].replace(",", ".")
                try:
                    pid = float(pid_str)
                except ValueError:
                    continue
                qstart = int(cols[idx["qstart"]])
                qend = int(cols[idx["qend"]])
                # Transform coordinates for flipped reads
                if flipped:
                    # original BLAST coords are 1-based on the unflipped read
                    qstart_new = read_len - qend + 1
                    qend_new = read_len - qstart + 1
                    qstart, qend = qstart_new, qend_new

                start = min(qstart, qend) - 1  
                end = max(qstart, qend)        
                key = (start, end)
                best = hits.get(key)
                if best is None or pid > best["pid"]:
                    hits[key] = {"pid": pid, "region": region}
        else:
            # headerless BLAST outfmt 6 
            # "read-id", "region", "percent-identity", "region-length", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "e-value", "bitscore"
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # whitespace split 
                cols = line.split()
                if len(cols) < 12:
                    continue
                qseqid = cols[0].strip()
                if qseqid != target_read_id:
                    continue
                sseqid  = cols[1]
                pident  = cols[2]
                qstart  = cols[7]
                qend    = cols[8]
                region = sseqid.strip().rstrip(".")
                pid = float(pident.replace(",", "."))
                qstart = int(qstart)
                qend   = int(qend)
                # If flipped, transform the coordinates on the flipped FASTA
                if flipped:
                    # if coords exceed read_len, something is wrong; skip
                    if qstart > read_len or qend > read_len:
                        continue
                    qstart_new = read_len - qend + 1
                    qend_new   = read_len - qstart + 1
                    qstart, qend = qstart_new, qend_new
                start = min(qstart, qend) - 1
                end   = max(qstart, qend)
                key = (start, end)
                best = hits.get(key)
                if best is None or pid > best["pid"]:
                    hits[key] = {"pid": pid, "region": region}
    if not hits:
        raise ValueError(f"No BLAST hits found for read {target_read_id} in {blast_path}")
    with open(bed_path, "w") as out:
        for (start, end) in sorted(hits.keys()):
            out.write(f"{target_read_id}\t{start}\t{end}\t{hits[(start,end)]['region']}\n")


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
                print(f"Standardfehlerausgabe: {result.stderr}")
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


def main():
    ap = argparse.ArgumentParser(
        description="""
        DUCKS4 - ID2bam2meth 
        
        Please choose for the reference either mode A or B:
        A) custom reference to be created from a read:
          --id_ref read-id --bam_ref read.bam --blast_ref blast.tsv/txt
        B) provide existing reference:
          --ref ref.fasta

        Alignment:
          --bam reads.bam
          [--id] read-id.txt   (if --id is omitted, whole bam is aligned. Here you can submit a chosen subset of reads to be filtered from the bam and aligned.)
        Methyl:
          --methyl
            - if custom reference: stats over BED from BLAST
            - if existing reference: require --region chr:start-end (creates 1-region BED) OR --regions_bed if more regions should be considered
        """
    )

    # Reference
      # ModeA
    ap.add_argument("--id_ref", required=False, help="Read-id used as custom reference (Mode A).")
    ap.add_argument("--bam_ref", required=False, help="BAM containing the reference-read (Mode A).")
    ap.add_argument("--blast_ref", required=False, help="BLAST CSV/TXT for reference-read annotation (folder blast_results) (Mode A).")
      # ModeB
    ap.add_argument("--ref", required=False, help="Use an existing reference FASTA (Mode B).")

    # optional reads for alignment
    ap.add_argument("--bam", required=False, help="Input BAM to extract/align reads against the reference.")
    ap.add_argument("--txt", required=False, help="Optional read_id.txt to subset reads from --bam to be aligned to reference.")

    # methylation
    ap.add_argument("--methyl", required=False, action='store_true', help="Run modkit pileup + stats.")
    ap.add_argument("--region", required=False, help="Required for Mode B: --ref + --methyl. Format: chr:start-end")
    ap.add_argument("--regions_bed", required=False, help="Optional BED with multiple regions for modkit stats (works for Mode A and B). If set, --region is not needed.")

    # outputs
    ap.add_argument("--out_prefix", required=False, default=None, help="Prefix name (default: id_ref or ref basename).")
    ap.add_argument("--out_path", required=False, help="Output directory. Default: created next to BAM.")
    ap.add_argument("--threads", dest="threads", default="45", required=False, help="Threads (default 45).")

    args = ap.parse_args()

    print("       ")
    print("####   DUCKS4 - ID2bam2meth   ####")
    print("       ")

    modeA = bool(args.id_ref and args.bam_ref and args.blast_ref)
    modeB = bool(args.ref)

    if not modeA and not modeB:
        print("ERROR: Either choose a read to create a custom reference from it OR give a reference .fasta file.")
        print("       custom reference needs: --id_ref read-id --bam_ref read.bam --blast_ref blast.tsv/txt from DUCKS4-output (folder blast-results)")
        print("       provided reference needs: --ref ref.fasta")
        sys.exit(1)

    if modeA and modeB:
        print("ERROR: Please choose either a custom reference (--id_ref, --bam_ref, blast_ref) to be created from a read OR provide an existing reference.fasta (..ref), not both.")
        sys.exit(1)

    if args.txt and not args.bam:
        print("ERROR: --bam is required when --txt is provided as reads are being filtered from the bam file.")
        sys.exit(1)

    if args.methyl and modeB and (not args.region) and (not args.regions_bed):
      print("ERROR: For --ref + --methyl you must provide either --region chr:start-end OR --regions_bed file.bed for methylation calculation.")
      sys.exit(1)

    if args.bam is None and args.methyl:
        print("ERROR: --bam is required if you want to run --methyl.")
        sys.exit(1)

    if args.out_prefix:
        prefix_name = args.out_prefix
    else:
        if modeA:
            prefix_name = args.id_ref
        else:
            prefix_name = Path(args.ref).stem

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

        if modeA:
            dir_name = f"custref2bam_{prefix_name}"
        else:
            dir_name = f"ref2bam_{prefix_name}"

        outpath = base_dir / dir_name
        outpath.mkdir(parents=True, exist_ok=True)

    def minimap2(file, ref_fasta, out_dir):
        print("       ")
        print("Alignment with Minimap2:")
        print("       ")
        file_in = os.path.basename(file)
        ref_name = os.path.basename(ref_fasta).split('.')[0]
        sam_file = ''.join([file_in.split('.')[0], "_", ref_name, ".sam"])

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

    # reference
    ref_fasta = None
    bed_path = None

    if modeA:
        bam_ref = Path(args.bam_ref)
        blast_ref = Path(args.blast_ref)

        fasta_path = outpath / f"{prefix_name}.fa"
        bed_path = outpath / f"{prefix_name}.bed"

        print("Create custom reference from read:")
        print(f"  id_ref    : {args.id_ref}")
        print(f"  bam_ref   : {bam_ref}")
        print(f"  blast_ref : {blast_ref}")

        seq, is_reverse, read_len, is_unmapped = get_read_sequence(bam_ref, args.id_ref)

        flipped = False
        if is_unmapped:
            print(f"Read {args.id_ref} is unaligned therefore read-orientation can't be determined.")
        elif is_reverse:
            print(f"Read {args.id_ref} is aligned reverse to reference – flipping to forward orientation.")
            seq = reverse_complement(seq)
            flipped = True
        else:
            print(f"Read {args.id_ref} is aligned in forward orientation.")
        

        write_fasta(args.id_ref, seq, fasta_path)
        index_fasta_with_samtools(fasta_path)
        parse_blast_to_bed(blast_ref, args.id_ref, bed_path, read_len, flipped)

        ref_fasta = fasta_path

        print("Done creating custom reference.")
        print(f"FASTA : {fasta_path}")
        print(f"FAIDX : {fasta_path}.fai")
        print(f"BED   : {bed_path}")

    else:
        # Mode B
        ref_fasta = Path(args.ref).resolve()
        print("Using provided reference FASTA:")
        print(f"  REF : {ref_fasta}")

        # optional: create faidx if missing
        if not Path(str(ref_fasta) + ".fai").exists():
            print("No .fai found -> indexing reference with samtools faidx.")
            index_fasta_with_samtools(ref_fasta)

    # aligning reads to reference
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

        # cleanup fastq if it was generated from bam
        if map_input.endswith(".bam") and fastq.endswith(".fastq") and Path(fastq).exists():
            subprocess.call(["rm", fastq])

        print(f"Aligned BAM: {outpath / bam_map}")

    else:
        bam_map = None
        print("No alignment requested (no --bam).")

    # Methylation with modkit
    if args.methyl:
        if bam_map is None:
            print("ERROR: --methyl requested but no aligned BAM exists.")
            sys.exit(1)

        print("       ")
        print("Start Methylation-analysis with Modkit.")
        print("       ")

        bam_map_name = Path(bam_map).stem
        meth_path = outpath / f"methylation-analysis_{prefix_name}"
        meth_path.mkdir(parents=True, exist_ok=True)

        modkit_bed = f"{bam_map_name}_modkit-methyl.bed"
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

        # regions to summarize:
        if args.regions_bed:
          regions_bed = Path(args.regions_bed).resolve()
        else:
          if modeA:
              regions_bed = bed_path
          else:
              # create a .bed from args.region
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

        # create bedGraph + bed from modkit stats 
        output_bg = f"{bam_map_name}.bedgraph"
        output_bed = f"{bam_map_name}.bed"

        chrom_idx, start_idx, end_idx, value_idx = 0, 1, 2, 10

        with open(modout_all, 'r') as fin, open(meth_path / output_bg, 'w', newline='') as fout:
            reader = csv.reader(fin, delimiter='\t')
            writer = csv.writer(fout, delimiter='\t')
            fout.write("track type=bedGraph name='Percent_Methylation' description='Percent Methylation'\n")
            next(reader, None)
            for row in reader:
                if len(row) <= value_idx:
                    continue
                bedgraph_row = [row[chrom_idx], row[start_idx], row[end_idx], row[value_idx]]
                writer.writerow(bedgraph_row)

        with open(modout_all, 'r') as fin, open(meth_path / output_bed, 'w', newline='') as fout:
            reader = csv.reader(fin, delimiter='\t')
            writer = csv.writer(fout, delimiter='\t')
            fout.write("track type=bed name='Methylation_Labels' description='Percent Methylation as Name'\n")
            next(reader, None)
            for row in reader:
                if len(row) <= value_idx:
                    continue
                try:
                    formatted_percent = f"{float(row[value_idx]):.2f}"
                    bed_row = [row[chrom_idx], row[start_idx], row[end_idx], formatted_percent]
                    writer.writerow(bed_row)
                except ValueError:
                    continue

        subprocess.call(["chown", "-R", "777", str(meth_path)])

        print("Methylation outputs:")
        print(f"  {modout_all}")
        print(f"  {meth_path / output_bg}")
        print(f"  {meth_path / output_bed}")

    print("       ")
    print("Workflow finished. Thank you for using this pipeline. TL :-)")
    print("       ")


if __name__ == "__main__":
    main()
