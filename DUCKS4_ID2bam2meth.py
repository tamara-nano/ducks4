import os, subprocess, time, shutil, sys
from argparse import ArgumentParser
from pprint import pprint


parser = ArgumentParser(
    description="""
    DUCKS4 - ID2bam2meth

    This script allows bundling of selected read IDs from DUCKS4 output
    for methylation analysis and BAM extraction using reference genome T2T-chm13v2.0. Also an alternative reference can be used.
    """
)
parser.add_argument("--txt", dest="tinput", help="required, read_id.txt, Copy the read IDs you want to bundle from the analysis files into a .txt file and give it to me :)", required=True)
parser.add_argument("--bam", dest="binput", help="required, provide mapped & sorted .bam file (e.g. from DUCKS4-output) for reference T2T-chm13v2.0. If you want to use another ref. Please add flag --ref", required=True)
parser.add_argument("--ref", dest="ref", help="optional, custom reference .fasta. Default is T2T_chm13v2.0.", required=False)
parser.add_argument("--methyl", dest="methyl", help="optional, perform methylation calling on selected reads.", required=False, action='store_true')
parser.add_argument("--region", dest="region", help=" optional, enomic region (e.g. chr1:1-100). REQUIRED if --ref & --methyl are set. Default for T2T_chm13v2.0 ref (when no --ref is given) = chr4:193540172-193543634 (2 most distal RU + gene-body).", required=False)
parser.add_argument("--threads", dest="threads", help="optional. Set your amount of threads. Default is 45", default="45", required=False)

parser.set_defaults(methyl=False)


args = parser.parse_args()

if args.ref and args.methyl and not args.region:
    print("Error: --region is required when --ref and --methyl are given.")
    sys.exit(1)
    
txt_path = os.path.dirname(args.tinput)
txt_name = os.path.basename(args.tinput).split('.')[0]
txt_file = os.path.basename(args.tinput)

bam_path = os.path.dirname(args.binput)
bam_name = os.path.basename(args.binput).split('.')[0]
bam_file = os.path.basename(args.binput)

script_path = os.path.abspath(os.path.dirname( __file__ ))


def main():

    print("       ")
    print("####   DUCKS4 - IDtobamtometh   ####")
    print("       ")
    print("DUCKS4 - Filter and map your read-ID.txt files from a bam alignment, and optional methylation basecalling.")
    print("       ")
    print("See --help for more infos.")
    print("The whole workflow is based on T2T-chm13v2.0 but other reference can be used.")
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
            print(f"{result.stderr}")
        except subprocess.CalledProcessError as e:
          print(f"Error while converting file: {e}")
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
          print(f"Errer while converting the file: {e}")
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
      ref_path = os.path.join(script_path, "ressources","reference", '')
      sam_file = ''.join([file_in.split('.')[0], "_T2Tchm13v2.sam"])
      if args.threads:
        if args.ref:
          ref_name = os.path.basename(args.ref).split('.')[0]
          sam_file = ''.join([file_in.split('.')[0], "_", ref_name, ".sam"])
          subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", args.threads, "-Y", "-y", "-o", os.path.join(path, sam_file), args.ref, os.path.join(path, file_in)])
        else:
          sam_file = ''.join([file_in.split('.')[0], "_T2Tchm13v2.sam"])
          subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", args.threads, "-Y", "-y", "-o", os.path.join(path, sam_file), os.path.join(ref_path, "chm13v2.0.fa"), os.path.join(path, file_in)])
      else:
        if args.ref:
          ref_name = os.path.basename(args.ref).split('.')[0]
          sam_file = ''.join([file_in.split('.')[0], "_", ref_name, ".sam"])
          subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", "45", "-Y", "-y", "-o", os.path.join(path, sam_file), args.ref, os.path.join(path, file_in)])
        else:
          sam_file = ''.join([file_in.split('.')[0], "_T2Tchm13v2.sam"])
          subprocess.call(["minimap2", "-ax", "lr:hq", "--MD", "-L", "-t", "45", "-Y", "-y", "-o", os.path.join(path, sam_file), os.path.join(ref_path, "chm13v2.0.fa"), os.path.join(path, file_in)])
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
    
    ## Filter reads from bam-alignment
    
    print("Filter reads from .bam-alignment")
    
    txt_path = os.path.dirname(args.tinput)
    txt_name = os.path.basename(args.tinput).split('.')[0]
    txt_file = os.path.basename(args.tinput)

    bam_path = os.path.dirname(args.binput)
    bam_name = os.path.basename(args.binput).split('.')[0]
    bam_file = os.path.basename(args.binput)
    
    ID_bam = bam_name + '_' + txt_name + '.bam'
    ID_out = os.path.join(bam_path, ID_bam)
    try:
      with open(ID_out, 'w') as out_f:
        result = subprocess.run(['samtools', 'view', '-h', '-N', args.tinput, args.binput], check=True, text=True, stdout=out_f, stderr=subprocess.PIPE)
      print(f"The file was successfully converted to {ID_out}.")
      if result.stderr:
        print(f"{result.stderr}")
    except subprocess.CalledProcessError as e:
      print(f"Error with converting file: {e}")
      print(f"{e.stderr}")
      return None
    
    map_file = checkfile_fastq(ID_out)
    sam_map = minimap2(map_file)
    bam_map = samtools_bam(os.path.join(bam_path, sam_map))

    
    ## Methylation calculation with modkit
    
    if args.methyl:
      print("       ")    
      print("Start Methylation-analysis with Modkit.")
      print("       ")
      
      if args.ref:
        ref_name = os.path.basename(args.ref).split('.')[0]
        ref_path = os.path.dirname(args.ref)
        ref_file = os.path.basename(args.ref)
      else:
        ref_path = os.path.join(script_path, "ressources","reference")
        ref_path = str(ref_path)
      
      bam_map_name = os.path.basename(bam_map).split('.')[0]
      meth_path = os.path.join(bam_path, "methylation-analysis")
      os.mkdir(meth_path)
      modkit_bed = ''.join([bam_map_name, "_", "modkit-methyl.bed"])   
      if args.ref:
        subprocess.call(["modkit", "pileup", os.path.join(bam_path, bam_map), os.path.join(meth_path, modkit_bed), "--cpg", "--ref", args.ref])
      else:
        subprocess.call(["modkit", "pileup", os.path.join(bam_path, bam_map), os.path.join(meth_path, modkit_bed), "--cpg", "--ref", os.path.join(ref_path, "chm13v2.0.fa")])
      modkit_bedgz = ''.join([modkit_bed, ".gz"])
      subprocess.call(["bgzip", os.path.join(meth_path, modkit_bed)])
      subprocess.call(["tabix", os.path.join(meth_path, modkit_bedgz)])
      
      if args.ref:
        chrom, coords = args.region.split(":")
        start, end = map(int, coords.split("-"))
        bed = ''.join(["coordinates_methcalc.bed"])
        bed_file = os.path.join(meth_path, bed)
        with open(bed_file, "w") as f:
          f.write(f"{chrom}\t{start - 1}\t{end}\n")
      else:
        if args.region:
          chrom, coords = args.region.split(":")
          start, end = map(int, coords.split("-"))
          bed = ''.join(["coordinates_methcalc.bed"])
          bed_file = os.path.join(meth_path, bed)
          with open(bed_file, "w") as f:
            f.write(f"{chrom}\t{start - 1}\t{end}\n")
        else:
          bed = ''.join(["distal-RU_gene-body_methcalc.bed"])
          bed_file = os.path.join(meth_path, bed)
          with open(bed_file, "w") as f:
            f.write("chr4\t193540172\t193543634\n")
        
      modout_all = "modkit-STATS.tsv" 
      subprocess.call(["modkit", "stats", "--regions", bed_file, "-o", os.path.join(meth_path, modout_all), os.path.join(meth_path, modkit_bedgz)])

    
      subprocess.call(["chown", "-R", "777", meth_path])

    print("       ")
    print("####   DUCKS4 - ID2bam2meth  ####")
    print("       ")
    print("Workflow is finished. Thank you for using this pipeline. :-) TL.")
    print("       ")
    print("####   DUCKS4 - ID2bam2meth   ####")
    print("       ")


if __name__ == "__main__":
    main()            
