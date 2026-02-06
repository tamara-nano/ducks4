# DUCKS4

FSHD-analysis tool for Nanopore-Sequencing.

Facioscapulohumeral Muscular Dystrophy (FSHD) is an autosomal dominant form of muscular dystrophy caused by genetic or epigenetic changes within the D4Z4-repeat at the DUX4-gene, on chromosome 4q. Genetic analysis is challenging due to a nearly identical region on chromosome 10, multiple haplotypes, long and short repeat subtypes, and complex rearrangements such as translocations and duplications. So far, no single method detects all known causes of FSHD.

We have developed an integrated approach combining an optimised wet-lab protocol with an automated bioinformatics workflow, called DUCKS4. It enables read-level resolution of the D4Z4 array for FSHD1 repeat sizing, variant detection for FSHD2, and detection of methylation patterns. Using NCBI BLAST, it assigns reads to chromosomes and haplotypes, supporting robust filtering and analysis. With long-read Nanopore sequencing technology, our tool enables precise determination of D4Z4 array size, individual haplotype assignment, methylation profiling, and complex allele analysis. It also allows for the detection of mosaicism and structural variation like interchromosomal translocations, providing a comprehensive, single-method solution for FSHD analysis.

**Update:** 
DUCKS4 also calls directly the PAS-sequence to each read if pLAM is available and provides it in the haplotypes-resolved output csv-files. It tests if the PAS is intact (4qA), 10qA or differently disrupted. 


## Prerequisite

This project requires [Docker](https://docs.docker.com/get-docker/) to be installed on your system.
Please follow the official installation instructions for your operating system.

## Installation

Pull docker image: 

`docker pull ghcr.io/tamara-nano/fshd_ducks4:latest`

or 

Build image with dockerfile:

download repository and unzip it \
`cd /path/ducks4/`

`docker build -t ghcr.io/tamara-nano/fshd_ducks4 .  `

## Usage

For running the tool:

`docker run --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest --input /data/mysample.bam --methyl --variant`

For showing more infos:

`docker run -it --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest --help`

| **tags** | **Infos** |
|:-----------------------------|:-----------------------------------------|
| --input | provide input-file. Best start with your basecalled SUP bam or fastq/fastq.gz-file. |
| --methyl | optional, methylation calling with modkit, target region: chr4:193540172-193543634. |
| --variant | optional, mapping to HG38, variant calling with clair3, sniffles2 (sniffles also calls against the T2T-aligned bam-file), phasing with whatshap and variant annotation with SNPEff and SnpSift against ClinVar (vv20250729). |
| --threads | optional, set threads. |
| --version | show program's version number and exit. |

The output is saved in the folder where the original input file is located.

## DUCKS4 output

DUCKS4 gives following output:

-   FSHD_overview-statistics: read and repeat counts separated for chromosomes and haplotypes in a csv.file
-   detailed statistics to each haplotype: csv-files separated into haplotypes with more detailed infos of reads and the resolved repeat-composition: ? Infos include: read.id, RU-count S//L, status (complete, partial read), resolved repeat-sequence, PAS-sequence in pLAM, PAS-type (4qA, 10qA, disrupted), warning
-   sorted and indexed bam files of reads sorted into haplotypes and aligned to T2T-chm13v2.0
-   sorted and indexed bam file of D4Z4-only-reads sorted to chr4 or 10 and reads aligned to T2T-chm13v2.0
-   sorted and indexed alignment.bam of all reads aligned to T2T-chm13v2.0
-   coverage.txt: coverage-infos for the alignment.bam called via samtools coverage, Coverage is calculated in the region chr4:192667301-192902247 upstream of the D4Z4-array.
-   folder with methylation-statistics and bed-methyl files for the 4qA haplotype and chimeric reads called with modkit.
-   folder with variant-calling results from clair3 and sniffles2 and annotated vcf.files (SNPEff and SnpSift against Clinvar database)
-   if variant calling is performed: haplotagged, sorted and index HG38 bam file for variant analysis
-   folder with original blast-results in csv-files, whereas the blast-output is also sorted to haplotypes

We recommend using the tables alongside with viewing the bam-files in a genome viewer like IGV-browser.


# Anaylsis of individual read-subsets

The DUCKS4-results make it easy to directly select reads for individal subset for further alignment, optional methylation-calling and analysis. 
If further subsets of reads should be filtered and analyzed. a read-id.txt needs to be provided along the alignment .bam-file.
Furthermore a custom reference can be created from a read-id or a existing reference can be given.

**Note**: Subsetting reads are f.ex. necessary when two 4qA alleles are present and the methylation status should be called. The tool can't distinguish two 4qA alleles and need to be further manually curated.

With this script it is possible to define a custom reference from a chosen read via read-id (f.ex. choose a complete read from the 4qA output). 
Optionally provide a subset reads to align against (either read-ids via TXT and/or a BAM-file) and optionally call methylation.
For the custom reference the blast-results of this reads are annotated within a annotation.bed file and if --methyl is chosen the average methylation will be calculated for each entry within the annotation.bed file. 
Therefore it is possible to get a methylation gradient-profile over the D4Z4-array if used with a custom reference and annotation.bed. 
The results can then further be inspected in a genome viewer like the IGV-browser.

## Mode A: for creating a custom reference & calling methylation over the D4Z4-array within the custom reference
`docker run --rm -v "$(pwd)":/data ghcr.io/tamara-nano/fshd_ducks4:latest id2bam2meth \
  --id_ref xxxx-xxxx-xxx-xxxx \   # input read_id \
  --bam_ref /data/sample.bam \
  --blast_ref /data/DUCKS4_output/blast-results/...fshd-blast.txt or .../...reads_blast.csv \
  --bam /data/sample.bam \
  --txt /data/read-ids.txt \
  --methyl 

**Note**: \
From the read-id a reference FASTA, FAIDX and annotation.bed file from the blast-output is created. \
Please be aware that the id_ref needs to be present in the bam_ref and blast_ref inputs. \
Methylation is called over all regions from the annotation.bed if no other regions are given (f.ex. --regions_bed). \
Therefore a methlyation-gradient over all D4Z4-RUs can be called and will be provided as .bedgraph output and for convenience as .bed file with values as labels.

## Mode B: providing an existing reference
`docker run --rm -v "$(pwd)":/data ghcr.io/tamara-nano/fshd_ducks4:latest id2bam2meth \
  --ref /data/ref.fasta \
  --bam /data/sample.bam \
  --txt /data/read-ids.txt \ 
  --methyl \
  --region chr:start-end \ OR
  --regions_bed data/regions.bed

**Note**: \
Methylation is called over the regions provided either as --region (one region) or --region_bed (several regions possible).

## Creating the read-ID.txt: 
Simply copy the reads-IDs you want to subset and filter from the DUCKS4-output tables into a txt-file:

Format read-id.txt:\

read-id1\
read-id3\
read-id5\
...

## For showing more infos:

`docker run -it --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest id2bam2meth --help`

| **tags** | **Infos** |
|:-----------------------------|:-----------------------------------------|
| **Mode A - create custom ref** |
| --id_ref | optional, Input read_id from a read (best a complete read) to use as the custom reference (query name) (Mode A).
| --bam_ref | optional, Input BAM file where the reference-read is located within (Mode A).
| --blast_ref | optional, BLAST TSV/TXT file containing the reference-read data (Mode A).
| --out_prefix | optional, Prefix for naming the reference (default: read-ID) (Mode A).
| **Mode B - existing ref** |
| --ref | optional, Use an existing reference FASTA (Mode B). |
| **reads for alignment** |
| --bam | optional, Input BAM to extract/align reads against the reference.
| --txt | optional, Optional read_id.txt to subset reads from --bam to be aligned to reference.
| **Methylation** |
| --methyl | optional, Mthylation calculation for the reads. Methylation calculation happens over the whole repeat-array for the mean CpG value for each D4Z4.
| --region | optional, Required for Mode B: --ref + --methyl. Format: chr:start-end.
| --regions_bed | optional, Optional BED with multiple regions for modkit stats (works for Mode A and B). If set, --region is not needed.
| **Outputs** |
| --threads | optional, Set your amount of threads. Default is 45.
| --out_path | optional, Give output_path, default: path from --bam_ref.

## Output:

- Mode A: reference.fasta, reference.fasta.fai, annotation.bed (from blast-output) generated from read_id
- aligned reads.bam/subset-reads.bam to reference
- Methylation: alignedreads.bedgraph, alignedreads.bed, alignedreads.methylbed, modkit-stats.tsv




# Further analysis

Sometimes it will be necessary to further determine the sub-haplotype of the allele. Therefore a scheme was developed to make it easy to distinguish the haplotypes (Tab.1) (sub-HP-help-sheet.xlsx). With the bed-file “Haplotypes_identification_regions.bed” (found in the DUCKS4 folder where also the script is) the necessary regions and also all relevant SNPs within D4F104S1 and the pLAM region as well as the restriction sites for BinI and XapI within the proximal D4Z4-RUs are marked. Relevant is the SSLP repeat in CLUHP-4-201 gene and the first 3-5 repeat units with the restriction enzyme sites of the D4Z4-array which needs to be manually inspected. There are three types of RU: chr4 – B-X+, chr10 – B+X- and a mix type – B-X-: The “+, plus” means the sequence of the restriction site is correct: B: CCTAGG and X: AAATTCC, if a SNP is found there then the restriction site is disabled “-, minus“. The blast-workflow only distinguishes between chr4 and chr10 repeat units and doesn't detect hybrid D4Z4 (B-X-) RU. The inspection of the restriction sites of the RU is only necessary in the case for 4A166Ha/b/c and 4A166 as the scheme itself is not enough to distinguish between those. 4A166 is NOT permissive for FSHD while 4A166H is! To further distinguish those haplotypes the analysis of the restriction sites for BinI (B) and XapI (X) is necessary. 4A166H has following D4Z4 order: c10-c10-c4…, while 4A166 has mix-mix-c4….

--> Haplotypes_identification_key.xlsx

## Example Data

Example data can be found on Figshare: 10.6084/m9.figshare.29930690

This repository contains the sequencing data from the human reference genomes HG001, HG002 and HG003 from the whole genome sequencing and adaptive sampling runs with long-read sequencing with Nanopore (Oxford Nanopore Technologies, UK). The high molecular weight (HMW) DNA from the cell-cultures were extracted with Monarch HMW-DNA Extraction Kit for Tissue (NEB, US) and the library prepared with SQK-ULK114 Kit (Oxford Nanopore Technologies, UK). The data were basecalled with Dorado basecaller with methylation calling for 5mCG and 5hmCG in SUP mode. The sequencing data are mapped, indexed and sorted bam files aligned to the T2T-chm13v2.0 reference and further filtered for the D4Z4 locus on chromosome 4 (4q35) and the homologous region on chromosome 10 (10q26). 



## Publication

If using the workflow for a publication please cite:

<Löwenstern T., Madritsch M., Horner D., Brait N., Güleray Lafci N., Schachner A., Gerykova Bujalkova M., Kałużewski T., Szyld P., Hengstschläger M., Dremsek P., Laccone F. DUCKS4: A comprehensive workflow for Nanopore sequencing analysis of Facioscapulohumeral Muscular Dystrophy (FSHD). Human Genomics. Manuscript accepted.>





