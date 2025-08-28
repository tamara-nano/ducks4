import os, subprocess, time, shutil, sys
import pysam

script_path = os.path.abspath(os.path.dirname( __file__ ))

input_file = sys.argv[1]  
print(f"PAS-check with aligned bam-file {input_file}")
FSHD_path =sys.argv[2]


output_txt = "PAS.txt"
output_txt = os.path.join(FSHD_path, output_txt)

chrom = "chr4"
start_pos = 193543619 # 1-based
length = 6
q4A_seq = "ATTAAA"
q10A_seq = "ATCAAA"

start_0 = start_pos - 1
end_0 = start_0 + length

bamfile = pysam.AlignmentFile(input_file, "rb")

rows = []

for read in bamfile.fetch(chrom, start_0, end_0):
    if read.is_secondary or read.is_unmapped:
        continue
    
    aligned_dict = dict(read.get_aligned_pairs(matches_only=False))
    reversed_dict = {v: k for k, v in aligned_dict.items()}
    bases = []
    for i in range(length):
      ref_pos = start_0 + i
      qpos = reversed_dict.get(ref_pos)
      base = read.query_sequence[qpos] if qpos is not None else "-"
      bases.append(base)

    read_seq = ''.join(bases)

    # Comparison to target
    if read_seq == q4A_seq:
      status = "4qA_PAS"
    elif read_seq == q10A_seq:
      status = "10qA_PAS"
    else:
      status = "PAS_disrupted"

    rows.append((read.query_name, read_seq, status))


## q10

chrom2 = "chr10"
start_pos2 = 134726542 # 1-basiert
length = 6
q4A_seq = "ATTAAA"
q10A_seq = "ATCAAA"

start_02 = start_pos2 - 1
end_02 = start_02 + length

bamfile = pysam.AlignmentFile(input_file, "rb")


for read in bamfile.fetch(chrom2, start_02, end_02):
    if read.is_secondary or read.is_unmapped:
        continue
    
    aligned_dict = dict(read.get_aligned_pairs(matches_only=False))
    reversed_dict = {v: k for k, v in aligned_dict.items()}
    bases = []
    for i in range(length):
      ref_pos = start_02 + i
      qpos = reversed_dict.get(ref_pos)
      base = read.query_sequence[qpos] if qpos is not None else "-"
      bases.append(base)

    read_seq = ''.join(bases)

    # Comparison to target
    if read_seq == q4A_seq:
      status = "4qA_PAS"
    elif read_seq == q10A_seq:
      status = "10qA_PAS"
    else:
      status = "PAS_disrupted"

    rows.append((read.query_name, read_seq, status))

with open(output_txt, "w") as f:
    f.write("read.id\tPAS.seq\tPAS.type\n")
    for row in rows:
        f.write("\t".join(row) + "\n")

print(f"File saved as '{output_txt}'")
