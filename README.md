# DUCKS4

FSHD-analysis tool for Nanopore-Sequencing.

Facioscapulohumeral Muscular Dystrophy (FSHD) is an autosomal dominant form of muscular dystrophy caused by genetic or epigenetic changes within the D4Z4-repeat at the DUX4-gene, on chromosome 4q. Genetic analysis is challenging due to a nearly identical region on chromosome 10, multiple haplotypes, long and short repeat subtypes, and complex rearrangements such as translocations and duplications. So far, no single method detects all known causes of FSHD.

We have developed an integrated approach combining an optimised wet-lab protocol with an automated bioinformatics workflow, called DUCKS4. It enables read-level resolution of the D4Z4 array for FSHD1 repeat sizing, variant detection for FSHD2, and detection of methylation patterns. Using NCBI BLAST, it assigns reads to chromosomes and haplotypes, supporting robust filtering and analysis. With long-read Nanopore sequencing technology, our tool enables precise determination of D4Z4 array size, individual haplotype assignment, methylation profiling, and complex allele analysis. It also allows for the detection of mosaicism and structural variation like interchromosomal translocations, providing a comprehensive, single-method solution for FSHD analysis.


**v1.1.0:**

Major changes:

Updated ID2bam2meth workflow:
- Removed dependency on external BLAST annotation files.
- Automatic reference-read re-BLASTing.
- Automatic orientation detection.
- Automatic annotation curation.
- Added distal_unit methylation region in ID2bam2meth.
- Simplified methylation workflow.
- Integrated per-RU BinI/XapI restriction site classification (DUCKS4 BX check) into the ID2bam2meth workflow. (Experimental)

Patch:

Fix dplyr compatibility.

**v1.0.1:** 

DUCKS4 directly determines the PAS sequence for each read if pLAM is present and reports PAS type (4qA, 10qA or disrupted) in the haplotype-resolved output tables.

## Prerequisite

This project requires [Docker](https://docs.docker.com/get-docker/) to be installed on your system.
Please follow the official installation instructions for your operating system.

## Installation

Pull docker image: 

`docker pull ghcr.io/tamara-nano/fshd_ducks4:v1.1.0`

or 

Build image with dockerfile:

download repository and unzip it \
`cd /path/ducks4/`

`docker build -t ghcr.io/tamara-nano/fshd_ducks4:v1.1.0 .  `

## Usage

For running the tool:

`docker run --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:v1.1.0 --input /data/mysample.bam --methyl --variant`

For showing more infos:

`docker run -it --rm -v $(pwd):/data ghcr.io/tamara-nano/fshd_ducks4:v1.1.0 --help`

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

| Output | Description |
|:-------|:------------|
| `FSHD_overview-statistics.csv` | Read and repeat counts per chromosome and haplotype |
| `{haplotype}_all-reads.csv` | Per-haplotype detailed read table including read.id, RU count, S/L, status (complete/partial), resolved repeat sequence, PAS sequence, PAS type (4qA, 10qA, disrupted), and warnings |
| `{haplotype}_complete-reads.csv` | Per-haplotype detailed read table for complete reads.
| `{haplotype}_{all/complete}_T2T.bam` / `.bai` | Sorted and indexed BAM of haplotype-assigned reads aligned to T2T-chm13v2.0 |
| `D4Z4_reads_chr4/chr10_T2T.bam` / `.bai` | Sorted and indexed BAM of D4Z4-only reads separated by chromosome, aligned to T2T-chm13v2.0 |
| `alignment_T2T.bam` / `.bai` | Sorted and indexed BAM of all reads aligned to T2T-chm13v2.0 |
| `coverage.txt` | Samtools coverage output for the region chr4:192667301-192902247 upstream of the D4Z4 array |
| `methylation/` | Methylation statistics and bedMethyl files for 4qA haplotype and chimeric reads called with modkit |
| `variant-calling/` | Clair3 and Sniffles2 results with annotated VCF files (SnpEff + SnpSift against ClinVar); includes haplotagged sorted and indexed HG38 BAM if variant calling was performed |
| `blast-results/` | Original BLAST output CSV files sorted by haplotype |

We recommend using the tables alongside with viewing the bam-files in a genome viewer like IGV-browser.


---

## Analysis of individual read-subsets - ID2BAM2METH

The DUCKS4-results make it easy to directly select reads for individual subsets for further alignment, optional methylation-calling and analysis. 
If further subsets of reads should be filtered and analyzed, a read-id.txt needs to be provided along the alignment .bam-file.
Furthermore a custom reference can be created from a read-id or an existing reference can be given.

The intended the workflow is: 
> **DUCKS4** → manual review and allele-level subsetting of results → **ID2bam2meth** for improved allele-specific alignment and methylation analysis

**Note**: Subsetting reads is necessary for example when two 4qA alleles are present and the methylation status should be called separately. The tool cannot distinguish two 4qA alleles and requires manual curation in that case.

With this script it is possible to define a custom reference from a chosen read via read-id (e.g. choose a complete read from the 4qA output). 
Optionally provide a subset of reads to align against (either read-ids via TXT and/or a BAM-file) and optionally call methylation.
For the custom reference the blast-results of this read are annotated within an annotation.bed file and if --methyl is chosen the average methylation will be calculated for each entry within the annotation.bed file. 
The results can then be further inspected in a genome viewer like the IGV-browser.

### Mode A: for creating a custom reference
```
docker run -it --rm -v "$(pwd)":/data ghcr.io/tamara-nano/fshd_ducks4:v1.1.0 id2bam2meth \
  --id_ref READ_ID \
  --bam_ref /data/sample-containing_id_ref-read.bam \
  --bam /data/sample.bam \
  --txt /data/read-ids.txt \
  --methyl
```

**Note**: \
From the read-id a reference FASTA, FAIDX and annotation.bed file from the blast-output is created. \
Please be aware that the id_ref needs to be present in the bam_ref input. \
Methylation is called over all regions from the annotation.bed if no other regions are given (e.g. --regions_bed). \
Therefore a methylation-gradient over all D4Z4-RUs can be called and will be provided as .bedgraph output and for convenience as .bed file with values as labels.
Also an aggregated region as distal_unit containing the most distal repeat and pLAM region is defined and methylation called for it.
The BX restriction site check is automatically run after methylation analysis and results are written to `{out_prefix}/D4Z4_BX_check/`.

### Mode B: providing an existing reference
```
docker run -it --rm -v "$(pwd)":/data ghcr.io/tamara-nano/fshd_ducks4:v1.2.0 id2bam2meth \
  --ref /data/ref.fasta \
  --bam /data/sample.bam \
  --txt /data/read-ids.txt \
  --methyl \
...
  optional:
...
  --region chr:start-end
```
or
```
  --regions_bed /data/regions.bed
```

**Note**: \
If no custom regions are provided, ID2bam2meth automatically generates annotations from the supplied reference FASTA using BLAST. \
The reference FASTA is automatically copied into the output folder for use with the BX check and IGV visualization.

### Creating the read-ID.txt: 
Simply copy the read-IDs you want to subset and filter from the DUCKS4-output tables into a txt-file:

Format read-id.txt:

```
read-id1
read-id3
read-id5
...
```

### For showing more infos:

`docker run -it --rm -v "$(pwd)":/data ghcr.io/tamara-nano/fshd_ducks4:v1.1.0 id2bam2meth --help`

| **tags** | **Infos** |
|:-----------------------------|:-----------------------------------------|
| **Mode A - create custom ref** | |
| --id_ref | optional, Input read_id from a read (best a complete read) to use as the custom reference (Mode A). |
| --bam_ref | optional, Input BAM file where the reference-read is located (Mode A). |
| --out_prefix | optional, Prefix for naming the reference (default: read-ID) (Mode A). |
| **Mode B - existing ref** | |
| --ref | optional, Use an existing reference FASTA (Mode B). |
| **reads for alignment** | |
| --bam | optional, Input BAM to extract/align reads against the reference. |
| --txt | optional, Optional read_id.txt to subset reads from --bam to be aligned to reference. |
| **Methylation** | |
| --methyl | optional, Methylation calculation for the reads. Methylation calculation happens over the whole repeat-array for the mean CpG value for each D4Z4. |
| --region | optional, single region in chr:start-end format. |
| --regions_bed | optional, BED file with multiple regions. |
| --skip-bx | optional, skip the BX restriction site check. |
| **Outputs** | |
| --threads | optional, Set your amount of threads. Default is 45. |
| --out_path | optional, Give output path, default: path from --bam_ref. |

### Output:

- Mode A: reference.fasta, reference.fasta.fai, annotation.bed generated from read_id
- Mode A & B: The annotation is generated automatically using BLAST and curate_annotation_bed.py.
- aligned reads.bam / subset-reads.bam to reference
- Methylation: alignedreads.bedgraph, alignedreads.bed, alignedreads.methylbed, modkit-stats.tsv for the annotated regions if no region/regions_bed is provided.
- D4Z4_BX_check/: BinI/XapI restriction site classification per RU (see BX check section above)

<p align="center">
  <a href="https://github.com/user-attachments/assets/ccc61679-d2fb-4f4a-bb27-b6028b02339c">
    <img src="https://github.com/user-attachments/assets/ccc61679-d2fb-4f4a-bb27-b6028b02339c" width="900">
  </a>
  
### Curated methylation output

Methylation values are reported for:

- CLUHP4
- D4F104S1
- individual D4Z4 repeat units
- pLAM
- DUX4_end

Additionally, the following aggregate regions are reported:

- distal_unit: terminal RU plus pLAM

distal_unit is included because distal methylation is particularly relevant for FSHD diagnostics.



---

### BinI/XapI per-RU restriction site classification (BX check) [EXPERIMENTAL]

The BX check is automatically run as part of the methylation workflow in ID2bam2meth and classifies each D4Z4 repeat unit (RU) of an assigned allele based on the presence or absence of two restriction enzyme recognition sequences detected directly from read alignments:

| Site | Enzyme | Sequence | D4Z4 type |
|:-----|:-------|:---------|:----------|
| B+ | BinI/AvrII | CCTAGG | intact in Chr10_D4Z4 |
| X+ | XapI/ApoI | AAATTCC | intact in Chr4_D4Z4 |

Three D4Z4 unit types are distinguished:

| Type | B site | X site | Meaning |
|:-----|:-------|:-------|:--------|
| Chr4_D4Z4 | B- | X+ | XapI site intact, BinI site absent |
| Chr10_D4Z4 | B+ | X- | BinI site intact, XapI site absent |
| Hybrid_D4Z4 | B- | X- | neither site intact — hybrid unit |

The sites are detected with up to 1 mismatch () tolerance (fuzzy matching) to account for SNVs that disrupt restriction sites. Exact matches and fuzzy matches (mm=1) are reported separately in the BED output.

Per-RU consensus is built by majority vote across all reads covering each RU position. Confidence = fraction of reads agreeing with the majority call. RUs with fewer than 3 covering reads are flagged as `low_coverage`; RUs where fewer than 60% of reads agree are flagged as `ambiguous_consensus`.

### BX check output

The BX check output is written to `{out_prefix}/D4Z4_BX_check/` and contains:

| File | Description |
|:-----|:------------|
| `BX_per_RU.csv` | Per-read per-RU classification table |
| `BX_summary.csv` | Per-read summary of RU type counts |
| `BX_consensus.csv` | Per-RU consensus type, confidence, and read counts |
| `array_structure.csv` | Compact array structure string (e.g. `c10-c4-c4-c4-c4-c4-c4-c4`) |
| `BX_sites.bed` | RU-level BED colored by D4Z4 type for IGV (itemRgb) |
| `BX_sites.bedgraph` | Per-RU consensus confidence as bedGraph for IGV |
| `BX_restriction_sites.bed` | BinI/XapI site positions (includes fuzzy matching) within each RU for IGV navigation |

IGV colors: blue = Chr10_D4Z4 (B+/X-), red = Chr4_D4Z4 (B-/X+), purple = Hybrid_D4Z4 (B-/X-), grey = ambiguous/unclassified, green distal D4Z4-S.

<p align="center">
  <img src="https://github.com/user-attachments/assets/7112765f-2af9-4706-9b18-55241017411c" width="900">
</p>


## Further analysis

Sometimes it will be necessary to further determine the sub-haplotype of the allele. Therefore a scheme was developed to make it easy to distinguish the haplotypes (Tab.1) (sub-HP-help-sheet.xlsx). With the bed-file "Haplotypes_identification_regions.bed" (found in the DUCKS4 folder where also the script is) the necessary regions and all relevant SNPs within D4F104S1 and the pLAM region as well as the restriction sites for BinI and XapI within the proximal D4Z4-RUs are marked. Relevant is the SSLP repeat in the CLUHP4-201 gene and the first 3–5 repeat units with the restriction enzyme sites of the D4Z4-array which need to be manually inspected in complex cases.

→ Haplotypes_identification_key.xlsx

## Example Data

Example data can be found on Figshare: 10.6084/m9.figshare.29930690

This repository contains the sequencing data from the human reference genomes HG001, HG002 and HG003 from whole genome sequencing and adaptive sampling runs with long-read sequencing with Nanopore (Oxford Nanopore Technologies, UK). The high molecular weight (HMW) DNA from the cell cultures were extracted with Monarch HMW-DNA Extraction Kit for Tissue (NEB, US) and the library prepared with SQK-ULK114 Kit (Oxford Nanopore Technologies, UK). The data were basecalled with Dorado basecaller with methylation calling for 5mCG and 5hmCG in SUP mode. The sequencing data are mapped, indexed and sorted BAM files aligned to the T2T-chm13v2.0 reference and further filtered for the D4Z4 locus on chromosome 4 (4q35) and the homologous region on chromosome 10 (10q26).

## Publication

If using the workflow for a publication please cite:

<Löwenstern T, Madritsch S, Horner D, Brait N, Güleray Lafci N, Schachner A, Gerykova Bujalkova M, Kałużewski T, Szyld P, Hengstschläger M, Dremsek P, Laccone F. DUCKS4: a comprehensive workflow for Nanopore sequencing analysis of facioscapulohumeral muscular dystrophy (FSHD). Hum Genomics. 2026 Feb 6;20(1):48. doi: 10.1186/s40246-026-00921-2. PMID: 41652446; PMCID: PMC12977661.>
