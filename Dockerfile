FROM ubuntu:22.04

LABEL org.opencontainers.image.title="DUCKS4"
LABEL org.opencontainers.image.description="FSHD analysis workflow for Nanopore reads"
LABEL org.opencontainers.image.source="https://github.com/tamara-nano/ducks4"
LABEL org.opencontainers.image.version="2.1.0"


# SYSTEM UPDATE & BASICS

ENV DEBIAN_FRONTEND=noninteractive
ENV CONDA_DIR=/opt/conda
ENV PATH=/usr/bin:$PATH

RUN apt-get update && apt-get install -y \
    wget git curl bzip2 unzip nano gcc g++ make cmake \
    python3 python3-pip python3-dev \
    r-base \
    zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libncursesw5 \
    libhts-dev libssl-dev libxml2-dev \
    libcurl4-gnutls-dev \
    openjdk-21-jre-headless \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
 

# PYTHON PACKAGES


COPY ressources/install/requirements.txt /tmp/

RUN /usr/bin/python3 -m pip install pysam
RUN /usr/bin/python3 -m pip install pandas
RUN /usr/bin/python3 -m pip install biopython
    
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


# ADD CHANNELS AND CREATE CLAIR3 ENVIRONMENT


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


# ENTRY


WORKDIR /ducks4
COPY . /ducks4/

RUN /usr/bin/java -Xmx4g -jar /ducks4/ressources/tools/snpEff/snpEff.jar download -v hg38 || true

RUN chmod +x /ducks4/DUCKS4.py
RUN chmod +x /ducks4/DUCKS4_ID2bam2meth.py

ENTRYPOINT ["/usr/bin/python3", "/ducks4/DUCKS4.py"]

CMD ["/bin/bash"]


