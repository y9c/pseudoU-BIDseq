# x86_64
FROM debian:bullseye-slim

ENV PATH="/opt/pipeline/bin:/bin:/opt/micromamba/bin:$PATH"
ENV MAMBA_ROOT_PREFIX="/opt/micromamba"
ENV MAMBA_EXE="/bin/micromamba"

# install system dependencies
RUN DEBIAN_FRONTEND=noninteractive apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install tzdata apt-utils wget curl git bzip2 make xsltproc gcc g++ pkg-config zlib1g-dev libxml2-dev python3 python3-distutils default-jre && apt-get clean && rm -rf /var/lib/apt/lists/*
# install package by micromamba
RUN wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C /tmp bin/micromamba && cp /tmp/bin/micromamba $MAMBA_EXE && $MAMBA_EXE shell init -s bash -p $MAMBA_ROOT_PREFIX && eval "$($MAMBA_EXE shell hook --shell=posix)" && \
  micromamba install -n base -y -c conda-forge -c bioconda libzlib falco bowtie2 star samtools bedtools fastp && micromamba clean -afy
# install package by python-pip
RUN pip install --no-cache-dir snakemake==7.16 cutadapt==4.1 panoptes-ui pysam parasail
# install umicollapse
RUN git clone --quiet --depth 1 https://github.com/Daniel-Liu-c0deb0t/UMICollapse.git /tmp/UMICollapse &&  cp /tmp/UMICollapse/umicollapse.jar /bin/umicollapse.jar && rm -rf /tmp/UMICollapse && \
  mkdir -p /bin/lib && wget -q -P /bin/lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && wget -q -P /bin/lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
# install cpup
RUN git clone --quiet --depth 1 https://github.com/y9c/cpup.git /tmp/cpup && make -C /tmp/cpup/ -j && cp /tmp/cpup/cpup /bin/cpup && rm -rf /tmp/cpup
# clean up and reduce size
RUN apt-get purge -y wget git bzip2 make xsltproc gcc g++ pkg-config && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN echo "This is a analysis pipelie (v1.0) for BID-seq data" > /README.md
ADD ./bin /opt/pipeline/bin
COPY ./Snakefile /opt/pipeline/Snakefile
COPY ./config.yaml /opt/pipeline/config.yaml
COPY ./bidseq /opt/pipeline/bin/bidseq
