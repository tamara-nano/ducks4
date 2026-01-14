FROM ubuntu:22.04

LABEL org.opencontainers.image.title="DUCKS4"
LABEL org.opencontainers.image.description="FSHD analysis workflow for Nanopore reads"
LABEL org.opencontainers.image.source="https://github.com/tamara-nano/ducks4"
LABEL org.opencontainers.image.version="2.3.0"


## SYSTEM UPDATE & BASICS

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=/usr/bin:$PATH

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget git curl bzip2 unzip nano gcc g++ make cmake \
    python3 python3-pip python3-dev \
    r-base \
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libncursesw5 \
    libhts-dev libssl-dev libxml2-dev \
    libcurl4-gnutls-dev \
    openjdk-21-jre-headless \
    ca-certificates \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
 

## PYTHON PACKAGES

COPY ressources/install/requirements.txt /tmp/

RUN /usr/bin/python3 -m pip install pysam pandas biopython
    
# Miniconda 

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p $CONDA_DIR && \
    rm /tmp/miniconda.sh && \
    $CONDA_DIR/bin/conda clean -afy
    
# R PACKAGES

RUN /usr/bin/Rscript -e "options(repos='https://cloud.r-project.org'); \
  install.packages('remotes'); \
  remotes::install_version('dplyr', version='1.1.4'); \
  remotes::install_version('tidyr', version='1.3.1')"
RUN /usr/bin/Rscript -e "library(dplyr); library(tidyr); cat('R OK\n')"


## BIOINFORMATICS TOOLS

# minimap2
RUN git clone https://github.com/lh3/minimap2 /tmp/minimap2 && \
    make -C /tmp/minimap2 && \
    cp /tmp/minimap2/minimap2 /usr/local/bin/ && \
    rm -rf /tmp/minimap2

# htslib 
RUN wget -q https://github.com/samtools/htslib/releases/download/1.20/htslib-1.20.tar.bz2 && \
    tar -xjf htslib-1.20.tar.bz2 && \
    cd htslib-1.20 && ./configure --prefix=/usr/local && make -j && make install && \
    cd / && rm -rf htslib-1.20 htslib-1.20.tar.bz2

# samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && ./configure --prefix=/usr/local && make -j && make install && \
    cd / && rm -rf samtools-1.20 samtools-1.20.tar.bz2

# seqtk
RUN git clone https://github.com/lh3/seqtk /tmp/seqtk && \
    make -C /tmp/seqtk && \
    cp /tmp/seqtk/seqtk /usr/local/bin/ && \
    rm -rf /tmp/seqtk

# modkit
ARG MODKIT_VERSION=v0.5.0
ARG MODKIT_ASSET=modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz
ARG MODKIT_URL=https://github.com/nanoporetech/modkit/releases/download/${MODKIT_VERSION}/${MODKIT_ASSET}

RUN set -eu; \
    work=/tmp/modkit_install; mkdir -p "$work"; \
    curl -fL "$MODKIT_URL" -o "$work/modkit.tgz"; \
    # optional: peek inside to confirm path; remove after it works once
    tar -tzf "$work/modkit.tgz" | head -n 20; \
    tar -xzf "$work/modkit.tgz" -C "$work"; \
    # find the binary wherever the archive put it
    BIN="$(find "$work" -type f -name modkit -perm -u+x | head -n1)"; \
    if [ -z "$BIN" ]; then echo "ERROR: modkit binary not found"; ls -R "$work"; exit 1; fi; \
    install -m 0755 "$BIN" /usr/local/bin/modkit; \
    rm -rf "$work"; \
    modkit --version
    
ENV PATH="/usr/local/bin:$PATH"

# create clair3 environment

SHELL ["/bin/bash", "-c"]

RUN /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
 
RUN /opt/conda/bin/conda install -y -n base conda-libmamba-solver && /opt/conda/bin/conda config --set solver libmamba

RUN /opt/conda/bin/conda create -y -n clair3 -c bioconda -c conda-forge -c defaults \
      python=3.9.0 \
      clair3
RUN /opt/conda/bin/conda clean -afy

RUN /opt/conda/bin/conda run -n clair3 pip install whatshap --no-cache-dir -r /tmp/requirements.txt
RUN /opt/conda/bin/conda run -n clair3 pip install sniffles --no-cache-dir -r /tmp/requirements.txt

ENV PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin


# create workdir

WORKDIR /ducks4
COPY . /ducks4/

COPY ducks4 /usr/local/bin/ducks4
RUN sed -i 's/\r$//' /usr/local/bin/ducks4 \
 && chmod +x /usr/local/bin/ducks4 \
 && chmod +x /ducks4/DUCKS4.py /ducks4/DUCKS4_ID2bam2meth.py


# clair3 model

RUN set -eux; \
  clair_dir="/ducks4/ressources/tools/clair3"; \
  model_name="r1041_e82_400bps_sup_v500"; \
  model_url="https://cdn.oxfordnanoportal.com/software/analysis/models/clair3/${model_name}.tar.gz"; \
  model_dir="${clair_dir}/${model_name}"; \
  tmpdir="$(mktemp -d)"; \
  \
  mkdir -p "${clair_dir}"; \
  echo "Downloading Clair3 model: ${model_name}"; \
  curl -fL "${model_url}" -o "${tmpdir}/${model_name}.tar.gz"; \
  tar -xzf "${tmpdir}/${model_name}.tar.gz" -C "${tmpdir}"; \
  \
  # The archive usually contains a top-level folder (often 'models/<model_name>' or '<model_name>')
  # Find the actual extracted model folder and move it into the canonical location.
  extracted="$(find "${tmpdir}" -maxdepth 4 -type d -name "${model_name}" | head -n1)"; \
  test -n "${extracted}"; \
  rm -rf "${model_dir}"; \
  mkdir -p "$(dirname "${model_dir}")"; \
  mv "${extracted}" "${model_dir}"; \
  \
  printf "model_name=%s\nsource_url=%s\ndownloaded_utc=%s\n" \
    "${model_name}" "${model_url}" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    > "${model_dir}/model_source.txt"; \
  rm -rf "${tmpdir}"

# snpEff + SnpSift
 
ARG SNPEFF_URL=https://snpeff-public.s3.amazonaws.com/versions/snpEff_latest_core.zip

RUN set -eux; \
    mkdir -p /ducks4/ressources/tools/snpEff; \
    curl -fL --retry 7 --retry-all-errors --connect-timeout 20 --max-time 600 \
      "${SNPEFF_URL}" -o /tmp/snpeff.zip; \
    unzip -q /tmp/snpeff.zip -d /tmp; \
    cp -a /tmp/snpEff/* /ducks4/ressources/tools/snpEff/; \
    rm -rf /tmp/snpeff.zip /tmp/snpEff; \
    test -f /ducks4/ressources/tools/snpEff/snpEff.jar; \
    test -f /ducks4/ressources/tools/snpEff/SnpSift.jar
    
ENV PATH="/ducks4/ressources/tools/snpEff:${PATH}"

# download clinvar db

RUN set -eux; \
  outdir=/ducks4/ressources/tools/snpEff; \
  mkdir -p "$outdir"; \
  url="https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz"; \
  url_tbi="${url}.tbi"; \
  curl -fL "$url" -o "$outdir/clinvar_hg38.vcf.gz"; \
  curl -fL "$url_tbi" -o "$outdir/clinvar_hg38.vcf.gz.tbi"; \
  gunzip -t "$outdir/clinvar_hg38.vcf.gz"; \
  { \
    echo "source_url=$url"; \
    echo "source_url_tbi=$url_tbi"; \
    echo "original_remote_name=clinvar.vcf.gz"; \
    echo "downloaded_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"; \
  } > "$outdir/clinvar_source.txt"

# blast

ARG BLAST_VERSION=2.14.0
ARG BLAST_TAR=ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
ARG BLAST_URL=https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/${BLAST_TAR}
ARG BLAST_DIR=/ducks4/ressources/tools/ncbi-blast-${BLAST_VERSION}+

RUN set -eux; \
    mkdir -p /ducks4/ressources/tools; \
    wget -q "${BLAST_URL}" -O /tmp/blast.tgz; \
    tar -xzf /tmp/blast.tgz -C /ducks4/ressources/tools; \
    rm -f /tmp/blast.tgz; \
    # sanity check
    test -x "${BLAST_DIR}/bin/makeblastdb"; \
    "${BLAST_DIR}/bin/makeblastdb" -version

ENV PATH="/ducks4/ressources/tools/ncbi-blast-2.14.0+/bin:${PATH}"

# references

RUN set -eux; \
  refdir=/ducks4/ressources/reference; \
  mkdir -p "$refdir"; \
  \
  # GrChr38
  HG38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"; \
  hg38_gz=/tmp/GRCh38_no_alt_analysis_set.fna.gz; \
  curl -fL "$HG38_URL" -o "$hg38_gz"; \
  gzip -dc "$hg38_gz" > "$refdir/hg38_no_alt.fa"; \
  rm -f "$hg38_gz"; \
  samtools faidx "$refdir/hg38_no_alt.fa"; \
  { \
    echo "name: hg38_no_alt.fa"; \
    echo "source_url: $HG38_URL"; \
    echo "original_file: $(basename "$HG38_URL")"; \
    echo "downloaded_utc: $(date -u +%Y-%m-%d)"; \
    echo "sha256_fa: $(sha256sum "$refdir/hg38_no_alt.fa" | awk '{print $1}')"; \
  } > "$refdir/hg38_no_alt.source.txt"; \
  \
  # T2T-chm13v2.0
  HS1_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz"; \
  hs1_gz=/tmp/hs1.fa.gz; \
  curl -fL "$HS1_URL" -o "$hs1_gz"; \
  gzip -dc "$hs1_gz" > "$refdir/chm13v2.0.fa"; \
  rm -f "$hs1_gz"; \
  samtools faidx "$refdir/chm13v2.0.fa"; \
  { \
    echo "name: chm13v2.0.fa"; \
    echo "source_url: $HS1_URL"; \
    echo "original_file: $(basename "$HS1_URL")"; \
    echo "downloaded_utc: $(date -u +%Y-%m-%d)"; \
    echo "sha256_fa: $(sha256sum "$refdir/chm13v2.0.fa" | awk '{print $1}')"; \
  } > "$refdir/chm13v2.0.source.txt"

# ENTRY



ENTRYPOINT ["/usr/local/bin/ducks4"]

CMD ["help"]


