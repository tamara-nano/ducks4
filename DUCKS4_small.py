import os, subprocess, time, shutil, sys
from argparse import ArgumentParser
from pprint import pprint

VERSION = "2.1.0"

# Implement ArgParser
parser = ArgumentParser(
    description=f"""
    DUCKS4 - a comprehensive FSHD-analysis workflow for Nanopore-reads.
    Small version without variant-calling.

    DUCKS4 {VERSION}.
    """
)
parser.add_argument("--input", dest="finput", help="Choose .fasta, .fastq, .fastq.gz or .ubam/.bam input file. Best to simply start with your SUP-ubam file from basecalling.", required=True)
parser.add_argument("--ref_t2t", dest="reft2t", help="Please provide the T2T-chm13 reference (.fasta)", required=True)
parser.add_argument("--methyl", dest="methyl", help="Optional methylation calculation for all 4qA-reads with modkit", required=False, action='store_true')
parser.add_argument("--threads", dest="threads", help="Set your amount of threads. Default is 45", default="45", required=False)
parser.add_argument("--version", action="version", version=f"DUCKS4 {VERSION}")

parser.set_defaults(methyl=False)
parser.set_defaults(threads=False)



args = parser.parse_args()


path_sample = os.path.dirname(args.finput)
file_name = os.path.basename(args.finput).split('.')[0]

fastq_file = os.path.basename(args.finput)

script_path = os.path.abspath(os.path.dirname( __file__ ))
db_path = os.path.join(script_path, "ressources","blast_db", "FSHD-blast")
blast_path = os.path.join(script_path, "ressources","tools", "ncbi-blast-2.14.0+", "bin", "blastn")
ducks4_path = os.path.join(script_path, "ressources","tools", '') 
rscript = os.path.join(ducks4_path, "FSHD_analysis.R")

if not os.path.isfile(args.finput):
    print(f"ERROR: Input file {args.finput} not found.")
    sys.exit(1)

if args.threads:
    try:
        args.threads = int(args.threads)
    except ValueError:
        print("ERROR: --threads must be an integer.")
        sys.exit(1)

def validate_reference(ref_path: str, auto_index: bool = True) -> str:
    ref_path = os.path.abspath(ref_path)
    if not os.path.isfile(ref_path):
        sys.exit(f"[ERROR] Reference FASTA not found: {ref_path}")
    fai_path = ref_path + ".fai"
    if os.path.isfile(fai_path):
        return ref_path  # all good
    if not auto_index:
        sys.exit(f"[ERROR] Missing index: {fai_path}\n"
                  f"        Create it with: samtools faidx {ref_path}")
    try:
        print(f"[INFO] Building FASTA index (.fai) for: {ref_path}")
        subprocess.run(["samtools", "faidx", ref_path], check=True)
    except subprocess.CalledProcessError as e:
        sys.exit(f"[ERROR] samtools faidx failed for {ref_path} (exit {e.returncode}).")
    if not os.path.isfile(fai_path):
        sys.exit(f"[ERROR] Index creation did not produce: {fai_path}")
    return ref_path



def main():

    print("       ")
    print("####   DUCKS4   ####")
    print("       ")
    print("DUCKS4 - a comprehensive FSHD-analysis pipeline for long-reads.")
    print("DUCKS4 v2.0")
    print("       ")
    print("See --help for more infos.")
    print("The whole workflow is based on T2T-chm13v2.0.")
    print("Code written and adapted by TL :-)")
    print("       ")
    
    
    def checkfile_fastq(file_in):
      file_out = file_in.rsplit('.', 1)[0] + '.fastq'
      if file_in.endswith('bam'):
        print(".bam-file is being converted to .fastq-file with samtools.")
        print("       ")
        file_out = file_in.rsplit('.', 1)[0] + '.fastq'
        try:
          with open(file_out, 'w') as out_f:
            result = subprocess.run(['samtools', 'fastq', '-T', '*', file_in], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
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
            result = subprocess.run(['gunzip', '-c', file_in], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
          print(f"The file was successfully converted to {file_out}.")
          if result.stderr:
            print(f"Standardfehlerausgabe: {result.stderr}")
        except subprocess.CalledProcessError as e:
          print(f"Error while converting the file: {e}")
          print(f"{e.stderr}")
          return None
      elif file_in.endswith('fastq'):
        file_in = file_out
      elif file_in.endswith('fasta'):
        file_in = file_out
      elif file_in.endswith('fa'):
        file_in = file_out
      else:
        print(f"This file-format is not supported!")
        return None
      return file_out
    
    def checkfile_fasta(file_in):
      # Test, if file exists:
      file_out = file_in.rsplit('.', 1)[0] + '.fasta'
      if not os.path.isfile(file_in):
        raise FileNotFoundError(f"The file {file_in} was not found.")
      if file_in.endswith('bam'):
        print("       ")
        print(".bam-file is being converted to .fasta-file with samtools.")
        print("       ")
        file_out = file_in.rsplit('.', 1)[0] + '.fasta'
        try:
          with open(file_out, 'w') as out_f:
            result = subprocess.run(['samtools', 'fasta', file_in], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
          print(f"The file was successfully converted to {file_out}.")
          if result.stderr:
            print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
          print(f"Error while converting the file: {e}")
          print(f"{e.stderr}")
          return None
      elif file_in.endswith('gz'):
        print("       ")
        print(".fastq.gz-file is being converted to .fasta-file with seqtk.")
        print("       ")
        try:
          with open(file_out, 'w') as out_f:
            result = subprocess.run(['seqtk', 'seq', '-a', file_in], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
          print(f"The file was successfully converted to {file_out}.")
          if result.stderr:
            print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
          print(f"Error while converting the file: {e}")
          print(f"{e.stderr}")
          return None
      elif file_in.endswith('fastq'):
        print("       ")
        print(".fastq-file is being converted to .fasta-file with seqtk.")
        print("       ")
        try:
          with open(file_out, 'w') as out_f:
            result = subprocess.run(['seqtk', 'seq', '-a', file_in], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
          print(f"The file was successfully converted to {file_out}.")
          if result.stderr:
            print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
          print(f"Error while converting the file: {e}")
          print(f"{e.stderr}")
          return None
      elif file_in.endswith('fasta'):
        print("       ")
        print(".fasta-file recognized, continue with Blast.")
        print("       ")
        file_out = file_in
      elif file_in.endswith('fa'):
        print("       ")
        print(".fasta-file recognized, continue with Blast.")
        print("       ")
        file_out = file_in
      else:
        print(f"This file-format is not supported!")
        return None
      return file_out
    
    def minimap2(file):
      print("       ")
      print("Mapping with Minimap2:")
      print("       ")
      file_in = os.path.basename(file)
      path = os.path.dirname(file)
      #ref_path = os.path.join(script_path, "ressources","reference", '')
      sam_file = ''.join([file_in.split('.')[0], "_T2Tchm13v2.sam"])
      if args.threads:
        subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", args.threads, "-Y", "-y", "-o", os.path.join(path, sam_file), args.reft2t, os.path.join(path, file_in)])
      else:
        subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", "45", "-Y", "-y", "-o", os.path.join(path, sam_file), args.reft2t, os.path.join(path, file_in)])
      return sam_file 
    
    def samtools_bam(sam_input):
      print("       ")
      print("Sorting and indexing of .sam-file:")
      print("       ")
      sam_file = os.path.basename(sam_input)
      path = os.path.dirname(sam_input)
      bam_file = ''.join([sam_file.split('.')[0], ".bam"])
      if args.threads:
        thread = ''.join(["-@", args.threads])
        subprocess.call(["samtools", "view", thread, "-S", "-b", os.path.join(path, sam_file), "-o", os.path.join(path, bam_file)])
        subprocess.call(["samtools", "sort", "-m", "15G", "-o", os.path.join(path, bam_file), os.path.join(path, sam_file)])
      else:
        subprocess.call(["samtools", "view", "-@45", "-S", "-b", os.path.join(path, sam_file), "-o", os.path.join(path, bam_file)])
        subprocess.call(["samtools", "sort", "-@16", "-m", "15G", "-o", os.path.join(path, bam_file), os.path.join(path, sam_file)])
      subprocess.call(["samtools", "index", os.path.join(path, bam_file)])
      subprocess.call(["rm", os.path.join(path, sam_file)])
      return bam_file

      
    def check_ID(ID):
      ID_file = os.path.basename(ID)
      path = os.path.dirname(ID)
      print("       ")
      print("Filtering haplotypes and complete reads and mapping to T2T-chm13v2.0. for ID ", ID_file)
      print("       ")
      file_size = os.path.getsize((os.path.join(path, ID_file)))
      if file_size == 0:
        print("No IDs in file ", ID_file, ". Skipped & ID-file deleted.")
        subprocess.call(["rm", os.path.join(path, ID_file)])
        return None
      else:
        ID_sam = ''.join([file_name, "_", ID_file.split('.')[0], ".sam"]) 
        ID_bam = ''.join([file_name, "_", ID_file.split('.')[0], ".bam"])
        fsf = open(os.path.join(path, ID_bam), "w")
        subprocess.call(['samtools', 'view',  '-h', '-N', os.path.join(path, ID_file), os.path.join(path_sample, bam_map)], stdout = fsf)
        fsf.close()
        fastq = checkfile_fastq(os.path.join(path, ID_bam))
        subprocess.call(["rm", os.path.join(path, ID_bam)])
        sam = minimap2(fastq)
        bam = samtools_bam(os.path.join(path, sam))
        subprocess.call(["rm", fastq])
        subprocess.call(["mv", os.path.join(path, ID_file), os.path.join(path, "read-IDs")])
      return bam




    ### START WORKFLOW  ###
    
    
    # check ref for .fai file
    print("Check reference-fasta if .fai index is available, else .fai index will be created.")
    ref_t2t = validate_reference(args.reft2t)
    print(f"T2T reference is ok: {ref_t2t}")
    
    ## Blast reads against FSHD-database
    
    fastq_file = os.path.basename(args.finput)
    file_name = fastq_file.split('.')[0]
    blast_file = ''.join([fastq_file.split('.')[0], "_fshd-blast.txt"])
    fasta_file = checkfile_fasta(os.path.join(path_sample, fastq_file))
    FSHD_path = os.path.join(path_sample, ''.join(["FSHD-analysis_", file_name]), '')
    os.mkdir(FSHD_path)
    print("       ")
    print("Start Blast-analysis!")
    print("       ")
    subprocess.call([blast_path, "-db", db_path, "-query", fasta_file, "-num_threads", "32", "-out", os.path.join(FSHD_path, blast_file), "-outfmt", '6 qseqid sseqid pident slen length mismatch gapopen qstart qend sstart send evalue bitscore'])
    subprocess.call(["rm", fasta_file])
    
    ## map reads
    print("Map files to reference.")
    map_file = checkfile_fastq(os.path.join(path_sample, fastq_file))
    sam_map = minimap2(map_file)
    bam_map = samtools_bam(os.path.join(path_sample, sam_map))
    
    
    print("       ")
    print("Start FSHD-analysis!")
    print("       ")
    
    ##PAS-checK
    print("PAS-sequences are analyzed.")
    passcript = os.path.join(ducks4_path, "PAS_check.py")
    subprocess.run(["/usr/bin/python3", passcript, os.path.join(path_sample, bam_map), FSHD_path])
    pas_file = "PAS.txt"
    
    
    ## Provide blast_file to R-script  
    print("Running analysis.")
    subprocess.call(["/usr/bin/Rscript", rscript, os.path.join(FSHD_path, blast_file), os.path.join(FSHD_path, "PAS.txt"), FSHD_path])   
    
    
    
    # Filter out Haplotypes and map
    ID_path = os.path.join(FSHD_path, "read-IDs")
    os.mkdir(ID_path)
    q4A_all = check_ID(os.path.join(FSHD_path, "4qA_all-reads-ID.txt"))
    q4A_complete = check_ID(os.path.join(FSHD_path, "4qA_complete-reads-ID.txt"))
    chimeric = check_ID(os.path.join(FSHD_path, "chimeric-reads-ID.txt"))
    q4B_all = check_ID(os.path.join(FSHD_path, "4qB_all-reads-ID.txt"))
    q4B_complete = check_ID(os.path.join(FSHD_path, "4qB_complete-reads-ID.txt"))
    D4Z4_chr4 = check_ID(os.path.join(FSHD_path, "D4Z4-only_chr4-reads-ID.txt"))
    D4Z4_chr10 = check_ID(os.path.join(FSHD_path, "D4Z4-only_chr10-reads-ID.txt"))
    
    chr4_undefined = check_ID(os.path.join(FSHD_path, "chr4-undefined_all-reads-ID.txt"))
    chr10_all = check_ID(os.path.join(FSHD_path, "chr10_all-reads-ID.txt"))
    chr10_complete = check_ID(os.path.join(FSHD_path, "chr10_complete-reads-ID.txt"))
    
    
    ## Calculation of coverage in the region chr4:192667301-192902247
    print("Calculate coverage with samtools coverage for region chr4:192667301-192902247.")
    coverage = "coverage.txt"
    subprocess.call(["samtools", "coverage",  "-r", "chr4:192667301-192902247", os.path.join(path_sample, bam_map), "-o", os.path.join(FSHD_path, coverage)])
    
    
    ## Methylation calculation with modkit
    
    if args.methyl:
      print("       ")    
      print("Start Methylation-analysis with Modkit.")
      print("       ")
      q4A_bam_complete = next((f for f in os.listdir(FSHD_path) if "4qA_complete-reads-ID" in f and f.endswith(".bam")), None)
      q4A_bam_all = next((f for f in os.listdir(FSHD_path) if "4qA_all-reads-ID" in f and f.endswith(".bam")), None)
      chimeric = next((f for f in os.listdir(FSHD_path) if "chimeric" in f and f.endswith(".bam")), None)
      #ref_path = os.path.join(script_path, "ressources","reference")
      #ref_path = str(ref_path)
      meth_path = os.path.join(FSHD_path, "methylation-analysis")
      os.mkdir(meth_path)
      if q4A_bam_all:
        modkit_bed_all = "4qA-all-reads_modkit-methyl.bed"
        q4A_bam_all = str(q4A_bam_all)
        subprocess.call(["modkit", "pileup", os.path.join(FSHD_path, q4A_bam_all), os.path.join(meth_path, modkit_bed_all), "--cpg", "--ref", args.reft2t])
        modkit_bed_allgz = "4qA-all-reads_modkit-methyl.bed.gz"
        subprocess.call(["bgzip", os.path.join(meth_path, modkit_bed_all)])
        subprocess.call(["tabix", os.path.join(meth_path, modkit_bed_allgz)])
        
        bed = ''.join(["distal-RU_gene-body_methcalc.bed"])
        bed_file = os.path.join(meth_path, bed)
      
        with open(bed_file, "w") as f:
          f.write("chr4\t193540172\t193543634\n")
        
        modout_all = "4qA-all-reads_modkit-STATS.tsv" 
        subprocess.call(["modkit", "stats", "--regions", bed_file, "-o", os.path.join(meth_path, modout_all), os.path.join(meth_path, modkit_bed_allgz)])
    
      if q4A_bam_complete:
        modkit_bed_complete = "4qA-complete-reads_modkit-methyl.bed"
        q4A_bam_complete = str(q4A_bam_complete)
        subprocess.call(["modkit", "pileup", os.path.join(FSHD_path, q4A_bam_complete), os.path.join(meth_path, modkit_bed_complete), "--cpg", "--ref", args.reft2t])
        modkit_bed_completegz = "4qA-complete-reads_modkit-methyl.bed.gz"
        subprocess.call(["bgzip", os.path.join(meth_path, modkit_bed_complete)])
        subprocess.call(["tabix", os.path.join(meth_path, modkit_bed_completegz)])
        modout_complete = "4qA-complete-reads_modkit-STATS.tsv"
        subprocess.call(["modkit", "stats", "--regions", bed_file, "-o", os.path.join(meth_path, modout_complete), os.path.join(meth_path, modkit_bed_completegz)])
        
      if chimeric:
        modkit_bed_chimeric = "chimeric-reads_modkit-methyl.bed"
        chimeric = str(chimeric)
        subprocess.call(["modkit", "pileup", os.path.join(FSHD_path, chimeric), os.path.join(meth_path, modkit_bed_chimeric), "--cpg", "--ref", args.reft2t])
        modkit_bed_chimericgz = "chimeric-reads_modkit-methyl.bed.gz"
        subprocess.call(["bgzip", os.path.join(meth_path, modkit_bed_chimeric)])
        subprocess.call(["tabix", os.path.join(meth_path, modkit_bed_chimericgz)])
        bed = ''.join(["distal-RU_gene-body_methcalc.bed"])
        bed_file = os.path.join(meth_path, bed)
        with open(bed_file, "w") as f:
          f.write("chr4\t193540172\t193543634\n")
        modout_chimeric = "chimeric-reads_modkit-STATS.tsv"
        subprocess.call(["modkit", "stats", "--regions", bed_file, "-o", os.path.join(meth_path, modout_chimeric), os.path.join(meth_path, modkit_bed_chimericgz)])

      else:
        print("No 4qA-Haplotypes or chimeric-reads found, therefore no Methylation-analysis.")


    
    subprocess.call(["chown", "-R", "777", FSHD_path])

    print("       ")
    print("####   DUCKS4   ####")
    print("       ")
    print("Workflow is finished. Thank you for using this pipeline. :-) TL.")
    print("       ")
    print("####   DUCKS4   ####")
    print("       ")


if __name__ == "__main__":
               
    if "--version" in sys.argv:
        print(f"DUCKS4 version {VERSION}")
        sys.exit(0)
        
    main() 
