
import os, subprocess, time, shutil, sys
import csv, re

script_path = os.path.abspath(os.path.dirname( __file__ ))


def main():

    input_file = sys.argv[1]  
    print(f"bam file used for variant calling with clair3 & sniffles2: {input_file}")
    FSHD_path = sys.argv[2]
    ref_path = sys.argv[3]
    reference = os.path.join(ref_path, "Homo_sapiens_GRCh38_no_alt.fasta")
    bam_t2t = sys.argv[4]
    if len(sys.argv) > 5:
      threads = sys.argv[5]
    else:
      threads = "4"
    
    bam_file = os.path.basename(input_file)
    path_sample = os.path.dirname(input_file)
    variant_path = os.path.join(FSHD_path, "variant-calling")
    os.mkdir(variant_path)
    

    print("Running Clair3 for Nanopore Kit14_400bps - small variant calling.")
    
    file_name = bam_file.split('.')[0]
    input_bam = ''.join(["--bam_fn=", os.path.join(path_sample, bam_file)])
    ref_file = ''.join(["--ref_fn=", reference])
    platform = ''.join(["--platform=","ont",])
    clair3_path = os.path.join(variant_path, 'clair3', '')
    os.mkdir(clair3_path)
    sample_name = ''.join(["--sample_name=", file_name])
    output = ''.join(["--output=", clair3_path])
    phase_path = os.path.join(clair3_path, 'tmp', 'phase_output', 'phase_vcf', '')
    model_path = ''.join(["--model_path=", os.path.join(script_path, "clair3/r1041_e82_400bps_sup_v500/")])
    phasing = "--enable_phasing"
    whatshap = "--whatshap=/opt/conda/envs/clair3/bin/whatshap"
    cthreads = ''.join(["--threads=", threads])

    subprocess.call(["/opt/conda/envs/clair3/bin/run_clair3.sh", input_bam, ref_file, cthreads, model_path, platform, output, sample_name, phasing, whatshap])

    zygo_vcf = ''.join(["phased_merge_output", ".vcf"])
    
    print("Haplotagging of HG38-bam file with whatshap.")    
    haplotag_bam = ''.join([file_name, "_haplotagged", ".bam"])
    ps = ''.join([file_name, "_haploblocks", ".gtf"])
    ps_path = ''.join(["--gtf=", os.path.join(variant_path, ps)])
    subprocess.call(["/opt/conda/envs/clair3/bin/whatshap", "haplotag", "-o",  os.path.join(variant_path, haplotag_bam), "--reference", reference, os.path.join(clair3_path, ''.join([zygo_vcf, ".gz"])), os.path.join(path_sample, bam_file), "--output-threads=20", "--ignore-read-groups", "--output-haplotag-list",  os.path.join(variant_path, "haplotag-list.tsv")])
    subprocess.call(["samtools", "index", os.path.join(variant_path, haplotag_bam)])
    subprocess.call(["/opt/conda/envs/clair3/bin/whatshap", "stats", ps_path, os.path.join(clair3_path, ''.join([zygo_vcf, ".gz"]))])
        
         
    # SNIFFLES2 - Structural Variant calling 

    print("Running Sniffles2 - Structural Variant Caller for HG38.")
    sniffles_file = ''.join([haplotag_bam.split('.')[0], "_sniffles2_phased", ".vcf"])
    subprocess.call(["/opt/conda/envs/clair3/bin/sniffles", "--input", os.path.join(variant_path, haplotag_bam), "--vcf", os.path.join(variant_path, sniffles_file), "--phase", "--output-rnames"])
    
    print("Running Sniffles2 - Structural Variant Caller for T2T.")
    bam_t2t_name = os.path.basename(bam_t2t)
    sniffles_file2 = ''.join([bam_t2t_name, "_sniffles2", ".vcf"])
    subprocess.call(["sniffles", "--input", os.path.join(path_sample, bam_t2t_name), "--vcf", os.path.join(path_sample, sniffles_file2), "--output-rnames"])
 
    # Variant annotation
    
    print("Variant-annotation to ClinVar-database with SnpSift and effect prediction with SnpEff.")
    db_csv = os.path.join(script_path, "Variants_FSHD_DB_T2T.csv")  
    snv_vcf = os.path.join(clair3_path, "phased_merge_output.vcf.gz")
    sv_vcf = os.path.join(variant_path, sniffles_file)

    # Annotation with snpEFF
    print("Annotation of SNVs with SNPeFF.")
    vcf_phase_annot_file = ''.join(["phased_merge_HG38_SnpEff", ".vcf"])
    fsf = open(os.path.join(clair3_path, vcf_phase_annot_file), "w")
    subprocess.call(["java", "-jar", os.path.join(script_path, "snpEff/snpEff.jar"), "hg38", "-noStats","-canon", os.path.join(clair3_path, "phased_merge_output.vcf.gz")], stdout = fsf)
    fsf.close()
    vcf_phase_annot_filegz = ''.join(["phased_merge_HG38_SnpEff", ".vcf.gz"])
    fsf = open(os.path.join(clair3_path, vcf_phase_annot_filegz), "wb")
    subprocess.call(["bgzip", "-c", os.path.join(clair3_path, vcf_phase_annot_file)], stdout = fsf)
    fsf.close()
    subprocess.call(["tabix", "-p", "vcf", os.path.join(clair3_path, vcf_phase_annot_filegz)])
    
    # Annotation with snpSift and Clinvar
    print("Annotation of SNVs with SNPSift and Clinvar-db.")
    # Annotate with ClinVar
    vcf_phase_annot_file2 = "phased_merge_HG38_snpeff-clinvar.vcf"
    with open(os.path.join(clair3_path, vcf_phase_annot_file2), "w") as fsf:
        subprocess.call([
            "java", "-Xmx1g", "-jar",
            os.path.join(script_path, "snpEff/SnpSift.jar"),
            "annotate", "-v",
            os.path.join(script_path, "snpEff/clinvar_20250729.vcf.gz"),
            os.path.join(clair3_path, vcf_phase_annot_filegz)
        ], stdout=fsf)
    
    # Compress and index annotated VCF
    vcf_phase_annot_filegz2 = "phased_merge_HG38_Clinvar.vcf.gz"
    with open(os.path.join(clair3_path, vcf_phase_annot_filegz2), "wb") as fsf:
        subprocess.call([
            "bgzip", "-c",
            os.path.join(clair3_path, vcf_phase_annot_file2)
        ], stdout=fsf)
    
    subprocess.call([
        "tabix", "-p", "vcf",
        os.path.join(clair3_path, vcf_phase_annot_filegz2)
    ])
    
    # Filter for FSHD-relevant SNVs
    print("Filtering for FSHD-relevant SNVs...")

    fshd_vcf = os.path.join(clair3_path, "fshd_relevant.vcf")
    fshd_filter = (
        "((ANN[*].GENE = 'DUX4') | "
        "(ANN[*].GENE = 'SMCHD1') | "
        "(ANN[*].GENE = 'LRIF1') | "
        "(ANN[*].GENE = 'DNMT3B') | "
        "(ANN[*].GENE = 'TRIM43') | "
        "(ANN[*].GENE = 'CAPN3') | "
        "(ANN[*].GENE = 'VCP') | "
        "(CLNDN =~ 'Facioscapulohumeral'))")

    with open(fshd_vcf, "w") as fsf:
        subprocess.run([
            "java", "-jar", os.path.join(script_path, "snpEff/SnpSift.jar"),
            "filter", fshd_filter,
            os.path.join(clair3_path, vcf_phase_annot_filegz2)
        ], stdout=fsf, check=True)
    
    # Convert filtered VCF to TSV
    print("Creating TSV table from FSHD-relevant SNVs...")
    fshd_tsv = os.path.join(clair3_path, "fshd_relevant.tsv")
    


    
    def parse_ann_field(ann_field):
        """Parse ANN= field from snpEff annotation."""
        ann_entries = []
        for ann in ann_field.split(","):
            parts = ann.split("|")
            if len(parts) > 4:
                ann_entries.append({
                    "effect": parts[1],
                    "impact": parts[2],
                    "gene": parts[3],
                    "gene_id": parts[4]
                })
        return ann_entries
    
    def get_info_value(info_str, key):
        match = re.search(rf"{re.escape(key)}=([^;]+)", info_str)
        return match.group(1) if match else ""
    
    with open(fshd_vcf) as vcf, open(fshd_tsv, "w", newline="") as out:
        writer = csv.writer(out, delimiter="\t")
        writer.writerow(["CHROM", "POS", "REF", "ALT", "GENE", "EFFECT", "IMPACT", "CLNSIG", "CLNDN"])
    
        for line in vcf:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom, pos, _id, ref, alt, _qual, _filter, info = fields[:8]
    
            ann_field = get_info_value(info, "ANN")
            clnsig = get_info_value(info, "CLNSIG")
            clndbn = get_info_value(info, "CLNDN")
    
            ann_entries = parse_ann_field(ann_field)
            for ann in ann_entries:
                writer.writerow([
                    chrom, pos, ref, alt,
                    ann["gene"], ann["effect"], ann["impact"],
                    clnsig, clndbn
                ])
    
    print(f"[INFO] FSHD-relevant VCF saved to: {fshd_vcf}")
    print(f"[INFO] FSHD-relevant TSV saved to: {fshd_tsv}")




if __name__ == "__main__":
    main()
    
