FROM ubuntu:22.04


# metadata 

LABEL org.opencontainers.image.title="DUCKS4_small"
LABEL org.opencontainers.image.description="FSHD analysis workflow for Nanopore reads"
LABEL org.opencontainers.image.source="https://github.com/tamara-nano/ducks4"
LABEL org.opencontainers.image.version="2.1.0"


# SYSTEM UPDATE & BASICS

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH=/usr/local/bin:/usr/bin:/bin

RUN apt-get update && apt-get install -y --no-install-recommends \
    wget git curl bzip2 unzip nano gcc g++ make cmake \
    python3 python3-pip python3-dev \
    r-base \
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libncursesw5 \
    libhts-dev libssl-dev libxml2-dev \
    libcurl4-gnutls-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
 
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
    

# PYTHON PACKAGES

RUN /usr/bin/python3 -m pip install pysam
RUN /usr/bin/python3 -m pip install pandas
RUN /usr/bin/python3 -m pip install biopython


# R PACKAGES

RUN /usr/bin/Rscript -e "options(repos='https://cloud.r-project.org'); \
  install.packages('remotes'); \
  remotes::install_version('dplyr', version='1.1.4'); \
  remotes::install_version('tidyr', version='1.3.1')"
RUN /usr/bin/Rscript -e "library(dplyr); library(tidyr); cat('R OK\n')"

# BIOINFORMATICS TOOLS

# minimap2
RUN git clone https://github.com/lh3/minimap2 /tmp/minimap2 && \
    make -C /tmp/minimap2 && \
    cp /tmp/minimap2/minimap2 /usr/local/bin/ && \
    rm -rf /tmp/minimap2

# samtools
RUN wget -q https://github.com/samtools/samtools/releases/download/1.20/samtools-1.20.tar.bz2 && \
    tar -xjf samtools-1.20.tar.bz2 && \
    cd samtools-1.20 && ./configure --prefix=/usr/local && make -j && make install && \
    cd / && rm -rf samtools-1.20 samtools-1.20.tar.bz2

# htslib 
RUN wget -q https://github.com/samtools/htslib/releases/download/1.22.1/htslib-1.22.1.tar.bz2 && \
    tar -xjf htslib-1.22.1.tar.bz2 && \
    cd htslib-1.22.1 && ./configure --prefix=/usr/local && make -j && make install && \
    cd / && rm -rf htslib-1.22.1 htslib-1.22.1.tar.bz2

# seqtk
RUN git clone https://github.com/lh3/seqtk /tmp/seqtk && \
    make -C /tmp/seqtk && \
    cp /tmp/seqtk/seqtk /usr/local/bin/ && \
    rm -rf /tmp/seqtk

# modkit
ARG MODKIT_VERSION=v0.5.0
ARG MODKIT_ASSET=modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz
ARG MODKIT_URL=https://github.com/nanoporetech/modkit/releases/download/${MODKIT_VERSION}/${MODKIT_ASSET}

RUN set -eux; \
  work=/tmp/modkit_install; mkdir -p "$work"; \
  curl -fL "$MODKIT_URL" -o "$work/modkit.tgz"; \
  # test archive (no pipe to avoid SIGPIPE with pipefail)
  tar -tzf "$work/modkit.tgz" >/dev/null; \
  tar -xzf "$work/modkit.tgz" -C "$work"; \
  BIN="$(find "$work" -type f -name modkit -perm -u+x | head -n1)"; \
  test -n "$BIN"; \
  install -m 0755 "$BIN" /usr/local/bin/modkit; \
  rm -rf "$work"; \
  /usr/local/bin/modkit --version

ENV PATH="/usr/local/bin:$PATH"

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


# ENTRY


WORKDIR /ducks4
COPY . /ducks4/

RUN chmod +x /ducks4/DUCKS4_small.py
RUN chmod +x /ducks4/DUCKS4_ID2bam2meth.py

ENTRYPOINT ["/usr/bin/python3", "/ducks4/DUCKS4_small.py"]

CMD ["/bin/bash"]


