# DUCKS4
FSHD-analysis tool for ONT.

We have developed an integrated approach combining an optimized wet-lab protocol with an automated bioinformatics workflow, called DUCKS4. It determines pathogenic repeat changes for FSHD1, variant detection for FSHD2, and detection of methylation patterns by resolving the exact composition of the D4Z4-array on a read-based level. Followed by NCBI BLAST identification, the tool sorts the reads to the respective chromosomes and haplotype, and filters the reads accordingly, making the subsequent analysis easily accessible and systematic.

# Installation

docker pull ducks4

or 

git clone file.git

+ install requirements: minimap2, samtools, clair3, sniffles2, whatshap, pandas, pysam, R (with packages dplyr, tidyr)
  

# Usage

Show infos via:

python3 DUCKS4.py --help 

Best to simply start with your SUP-ubam file from basecalling. Or choose .fasta, .fastq, .fastq.gz or .ubam/.bam input file. The script will detect your file-type and converts it accordingly for what is needed.

python3 DUCKS4.py --input bam/fastq --methyl --variant --threads

--methyl     optional methylation basecalling with modkit
--variant    optional variant-calling and phasing with clair3, sniffles2 and whatshap
--threads    optional, set threads

# anaylsis of individual read-subsets

If further subsets of reads should be filtered and analyzed. a read-id.txt needs to be provided along the alignment .bam-file.

Simply copy the reads-IDs you want to subset and filter from the DUCKS4-output tables into a txt-file:

Format read-id.txt: 
read-id1
read-id3
read-id5
...

python3 DUCKS4_ID2bam2meth.py --help

python3 DUCKS4_ID2bam2meth.py --txt reads-d.txt --bam alignment.bam --methyl --regions

All parameters:
--txt        provide read-id.txt
--bam        provide alignment.bam
--ref        optional, provide own reference, if not the T2T-chm13v2.0 Ref from the DUCKs4-wf is used
--methyl     optional, methylation calling with modkit, target region: chr4:193540172-193543634
--region     optional, provide own region for methylation calling, format chr4:1-200
--threads    optional, set threads

The output is saved in the folder where the alignment.bam is located.

# Output

- FSHD_overview-statistics: read and repeat counts from the found haplotypes
- detailed statistics to each haplotype: sv-files separated into haplotypes with more detailed infos of reads and the resolved repeat-composition
- reads sorted into haplotypes and mapped to T2T-chm13v2.0
- D4Z4-only-reads sorted to chr4 or 10 and chr4 reads mapped to T2T-chm13v2.0
- alignment.bam of all reads aligned, sorted and indexed to T2T-chm13v2.0
- coverage.txt: coverage-infos for the alignment.bam called via samtools coverage, Coverage is calculated in the region chr4:192667301-192902247 upstream of the D4Z4-array.
- folder with methylation-statistics for the 4qA haplotype and chimeric reads called with modkit
- folder with variant-calling results
- folder with original blast-results, also sorted to haplotypes

# Further analysis

Sometimes it will be necessary to further determine the sub-haplotype of the allele. Therefore a scheme was developed to make it easy to distinguish the haplotypes (Tab.1) (sub-HP-help-sheet.xlsx). With the bed-file “Subhaplotypes-regions.bed” (found in the DUCKS4 folder where also the script is) the necessary regions and also all relevant SNPs within D4F104S1 and the pLAM region as well as the restriction sites for BinI and XapI within the proximal D4Z4-RUs are marked. Relevant is the SSLP repeat in CLUHP-4-201 gene and the first 3-5 repeat units with the restriction enzyme sites of the D4Z4-array which needs to be manually inspected. There are three types of RU: chr4 – B-X+, chr10 – B+X- and a mix type – B-X-: The “+, plus” means the sequence of the restriction site is correct: B: CCTAGG and X: AAATTCC, if a SNP is found there then the restriction site is disabled “-, minus“. The blast-workflow only distinguishes between chr4 and chr10 repeat units and doesn't detect hybrid D4Z4 (B-X-) RU.
The inspection of the restriction sites of the RU is only necessary in the case for 4A166Ha/b/c and 4A166 as the scheme itself is not enough to distinguish between those. 4A166 is NOT permissive for FSHD while 4A166H is! To further distinguish those haplotypes the analysis of the restriction sites for BinI (B) and XapI (X) is necessary. 4A166H has following D4Z4 order: c10-c10-c4…, while 4A166 has mix-mix-c4….

< Insert scheme here>

