# DUCKS4

FSHD-analysis tool for Nanopore-Sequencing.

Facioscapulohumeral Muscular Dystrophy (FSHD) is an autosomal dominant form of muscular dystrophy caused by genetic or epigenetic changes within the D4Z4-repeat at the DUX4-gene, on chromosome 4q. Genetic analysis is challenging due to a nearly identical region on chromosome 10, multiple haplotypes, long and short repeat subtypes, and complex rearrangements such as translocations and duplications. So far, no single method detects all known causes of FSHD.

We have developed an integrated approach combining an optimized wet-lab protocol with an automated bioinformatics workflow, called DUCKS4. It enables read-level resolution of the D4Z4 array for FSHD1 repeat sizing, variant detection for FSHD2, and detection of methylation patterns. Using NCBI BLAST, it assigns reads to chromosomes and haplotypes, supporting robust filtering and analysis. With long-read Nanopore sequencing technology, our tool enables precise determination of D4Z4 array size, individual haplotype assignment, methylation profiling, and complex allele analysis. It also allows for the detection of mosaicism and structural variation like interchromosomal translocations, providing a comprehensive, single-method solution for FSHD analysis.

**Update:** 
v2.1.0: DUCKS4 now also calls directly the PAS-sequence to each read if pLAM is available and provides it in the haplotypes-resolved output csv-files. It tests if the PAS is intact (4qA), 10qA or differently disrupted. 

## Installation

Pull docker image: 

docker pull fshd_ducks4

or 

Build image with dockerfile:

download repository and unzip it \
cd /path/ducks4/  \
docker build -t fshd_ducks4 .  

## Usage

For running the tool:

`docker run --rm -v $(pwd):/data fshd_ducks4 --input /data/mysample.bam --methyl --variant`

For showing more infos:

`docker run -it --rm -v $(pwd):/data fshd_ducks4 --help`

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

`docker run --rm -v "$(pwd)":/data --entrypoint python3 fshd-ducks4 /ducks4/DUCKS4_ID2bam2meth.py
  --txt /data/read-id.txt
  --bam /data/sample.bam
  --methyl
  --region chr4:193540172-193543634`

For showing more infos:

`docker run -it --rm -v $(pwd):/data fshd_ducks4 python3 /ducks4/DUCKS4_ID2bam2meth.py --help`

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

Sometimes it will be necessary to further determine the sub-haplotype of the allele. Therefore a scheme was developed to make it easy to distinguish the haplotypes (Tab.1) (sub-HP-help-sheet.xlsx). With the bed-file “Haplotypes_identification_regions.bed” (found in the DUCKS4 folder where also the script is) the necessary regions and also all relevant SNPs within D4F104S1 and the pLAM region as well as the restriction sites for BinI and XapI within the proximal D4Z4-RUs are marked. Relevant is the SSLP repeat in CLUHP-4-201 gene and the first 3-5 repeat units with the restriction enzyme sites of the D4Z4-array which needs to be manually inspected. There are three types of RU: chr4 – B-X+, chr10 – B+X- and a mix type – B-X-: The “+, plus” means the sequence of the restriction site is correct: B: CCTAGG and X: AAATTCC, if a SNP is found there then the restriction site is disabled “-, minus“. The blast-workflow only distinguishes between chr4 and chr10 repeat units and doesn't detect hybrid D4Z4 (B-X-) RU. The inspection of the restriction sites of the RU is only necessary in the case for 4A166Ha/b/c and 4A166 as the scheme itself is not enough to distinguish between those. 4A166 is NOT permissive for FSHD while 4A166H is! To further distinguish those haplotypes the analysis of the restriction sites for BinI (B) and XapI (X) is necessary. 4A166H has following D4Z4 order: c10-c10-c4…, while 4A166 has mix-mix-c4….

\< Insert scheme here\>

## Example Data

HG001

HG002

HG003


## Publication

If using the workflow for a publication please cite:

<Insert Paper-citation>

















