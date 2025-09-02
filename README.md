# DUCKS4

FSHD-analysis tool for Nanopore-Sequencing.

Facioscapulohumeral Muscular Dystrophy (FSHD) is an autosomal dominant form of muscular dystrophy caused by genetic or epigenetic changes within the D4Z4-repeat at the DUX4-gene, on chromosome 4q. Genetic analysis is challenging due to a nearly identical region on chromosome 10, multiple haplotypes, long and short repeat subtypes, and complex rearrangements such as translocations and duplications. So far, no single method detects all known causes of FSHD.

We have developed an integrated approach combining an optimized wet-lab protocol with an automated bioinformatics workflow, called DUCKS4. It enables read-level resolution of the D4Z4 array for FSHD1 repeat sizing, variant detection for FSHD2, and detection of methylation patterns. Using NCBI BLAST, it assigns reads to chromosomes and haplotypes, supporting robust filtering and analysis. With long-read Nanopore sequencing technology, our tool enables precise determination of D4Z4 array size, individual haplotype assignment, methylation profiling, and complex allele analysis. It also allows for the detection of mosaicism and structural variation like interchromosomal translocations, providing a comprehensive, single-method solution for FSHD analysis.

**Update:** 
v2.1.0: DUCKS4 now also calls directly the PAS-sequence to each read if pLAM is available and provides it in the haplotypes-resolved output csv-files. It tests if the PAS is intact (4qA), 10qA or differently disrupted. 

## Prerequisite

This project requires [Docker](https://docs.docker.com/get-docker/) to be installed on your system.
Please follow the official installation instructions for your operating system.

## Installation

Pull docker image: 

`docker pull ghcr.io/tamara-nano/fshd_ducks4:latest`

**Note:**

This repository only contains the scripts, it can't therefore not be build via dockerfile.

## Usage

For running the tool:

`docker run --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest --input /data/mysample.bam --methyl --variant`

For showing more infos:

`docker run -it --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest --help`

| **tags** | **Infos** |
|:-----------------------------|:-----------------------------------------|
| --bam | provide input-file. Best start with your basecalled SUP bam or fastq/fastq.gz-file. |
| --methyl | optional, methylation calling with modkit, target region: chr4:193540172-193543634. |
| --variant | optional, mapping to HG38, variant calling with clair3, sniffles2 (sniffles also calls against the T2T-aligned bam-file), phasing with whatshap and variant annotation with SNPEff and SnpSift against ClinVar (vv20250729). |
| --threads | optional, set threads. |

The output is saved in the folder where the original input file is located.

## Anaylsis of individual read-subsets

The DUCKS4-results makes it easy to directly select reads for individal subset for further mapping, optional methylation-calling and analysis. If further subsets of reads should be filtered and analyzed. a read-id.txt needs to be provided along the alignment .bam-file.

**Note**: Subsetting reads are f.ex. necessary when 2 4qA alleles are present and the methylation status should be called. The tool can't distinguish two 4qA alleles.

`docker run --rm -v "$(pwd)":/data --entrypoint python3 ghcr.io/tamara-nano/fshd_ducks4:latest /ducks4/DUCKS4_ID2bam2meth.py \
  --txt /data/read-id.txt \
  --bam /data/sample.bam \
  --methyl \
  --region chr4:193540172-193543634`

For showing more infos:

`docker run -it --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:latest python3 /ducks4/DUCKS4_ID2bam2meth.py --help`

| **tags** | **Infos** |
|:-----------------------------|:-----------------------------------------|
| --txt | required, provide read-id.txt: Copy the read IDs you want to bundle from the analysis files into a .txt file. |
| --bam | required, provide mapped & sorted .bam file for reference T2T-chm13v2.0 (e.g. from DUCKS4-output). If you want to use another ref. Please add flag --ref --ref optional, provide own reference, else the T2T-chm13v2.0 ref from the DUCKs4-wf is used. |
| --ref | optional, provide own reference, else the T2T-chm13v2.0 ref from the DUCKs4-wf is used. |
| --methyl | optional, methylation calling with modkit, target region: chr4:193540172-193543634 (2 most distal RU + gene-body with pLAM). |
| --region | optional, provide genomic region (e.g. chr1:1-100). Only REQUIRED if --ref & --methyl are set. Default for T2T_chm13v2.0 ref (when no --ref is given) = chr4:193540172-193543634. |
| --threads | optional, set threads. |

The output is saved in the folder where the alignment.bam is located.

**Attention:** If the flag --ref and --methyl are set, the flag --region must also be provided!

Creating the read-ID.txt Simply copy the reads-IDs you want to subset and filter from the DUCKS4-output tables into a txt-file:

Format read-id.txt:\

read-id1\
read-id3\
read-id5\
...

## DUCKS4 output

DUCKS4 gives following output:

-   FSHD_overview-statistics: read and repeat counts separated for chromosomes and haplotypes in a csv.file
-   detailed statistics to each haplotype: csv-files separated into haplotypes with more detailed infos of reads and the resolved repeat-composition: ? Infos include: read.id, RU-count S//L, status (complete, partial read), resolved repeat-sequence, PAS-sequence in pLAM, PAS-type (4qA, 10qA, disrupted), warning
-   sorted and indexed bam files of reads sorted into haplotypes and mapped to T2T-chm13v2.0
-   sorted and indexed bam file of D4Z4-only-reads sorted to chr4 or 10 and chr4 reads mapped to T2T-chm13v2.0
-   sorted and indexed alignment.bam of all reads aligned to T2T-chm13v2.0
-   coverage.txt: coverage-infos for the alignment.bam called via samtools coverage, Coverage is calculated in the region chr4:192667301-192902247 upstream of the D4Z4-array.
-   folder with methylation-statistics and bed-methyl files for the 4qA haplotype and chimeric reads called with modkit.
-   folder with variant-calling results from clair3 and sniffles2 and annotated vcf.files (SNPEff and SnpSift against Clinvar database)
-   if variant calling is performed: haplotagged, sorted and index HG38 bam file for variant analysis
-   folder with original blast-results in csv-files, whereas the blast-output is also sorted to haplotypes

## Further analysis

Sometimes it will be necessary to further determine the sub-haplotype of the alleles. Therefore we developed a haplotype identification key to make it easy to distinguish the haplotypes based on the results from the DUCKS4 workflow. With the bed-file “Haplotypes_identification_regions.bed” (found in the DUCKS4 folder where also the script is) the necessary regions are marked.
The table is to be read from left to right: the first 6 columns (title DUCKS4) are the results from the DUCKS4 detailed output for the haplotypes and the composition of the repeat-sequence. The columns 7 to 9 (title Haplotype results) shows the resulting haplotype, permissive haplotypes are marked red. The columns 10-11 (title manual check) are displaying the data relevant for the manual check. Relevant is the SSLP repeat in CLUHP-4-201 gene (SSLP: chr4:193430066-193430226, T2T-chm13v2.0) and the first 3-5 repeat units with the restriction enzyme sites of the D4Z4-array for manual inspection within the alignment bam-files of the haplotypes from the DUCKS4-output. With the restriction enzyme sites (BinI (B)/XapI (X)) from Southern Blotting three types of RU are discerned: chr4 – B-X+, chr10 – B+X- and a hybrid type – B-X-: The “+, plus” means the sequence of the restriction site is correct: BinI: CCTAGG and XapI: AAATTCC. If a SNP is found there then the restriction site is disrupted “-, minus“. The DUCKS4 blast-workflow only distinguishes between chr4 and chr10 repeat units, we experienced so far that the hybrid repeat units are called as chr4 types with DUCKS4 (but we haven't validated it). The hybrid type is only relevant for determining the non-permissive 4qA subtype 4A166 as its proximal repeat units start with hybrid RU (D4Z4-order: hyb – hyb – mix hyb/c4) but also the permissive 4A166H should be verified as its proximal repeat units start with chr10 RU (D4Z4-order: c10 – c10 – mix c4/c10). Please refer to Giardina et al. 2024 (details in the Supplementary Infos), as well as Lemmers et al. 2010a and 2022a for further information to the various haplotypes and hybrid alleles.

\< Insert scheme here\>

## Example Data

Example sequence data can be found on Figshare:
[10.6084/m9.figshare.29930690](https://figshare.com/projects/DUCKS4-FSHD/260306)

This repository contains the sequencing data from the human reference genomes HG001, HG002 and HG003 from the whole genome sequencing and adaptive sampling runs with long-read sequencing with Nanopore (Oxford Nanopore Technologies, UK). The high molecular weight (HMW) DNA from the cell-cultures were extracted with Monarch HMW-DNA Extraction Kit for Tissue (NEB, US) and the library prepared with SQK-ULK114 Kit (Oxford Nanopore Technologies, UK). The data were basecalled with Dorado basecaller with methylation calling for 5mCG and 5hmCG in SUP mode. The sequencing data are mapped, indexed and sorted bam files aligned to the T2T-chm13v2.0 reference and further filtered for the D4Z4 locus on chromosome 4 (4q35) and the homologous region on chromosome 10 (10q26). 



## Publication

If using the workflow for a publication please cite:

<Löwenstern T., Madritsch M., Horner D., Brait N., Güleray Lafci N., Schachner A., Gerykova Bujalkova M., Kałużewski T., Szyld P., Hengstschläger M., Dremsek P., Laccone F. DUCKS4: A comprehensive workflow for Nanopore sequencing analysis of Facioscapulohumeral Muscular Dystrophy (FSHD). Manuscript in preparation.>

























