FROM ubuntu:20.04
RUN echo "This is a analysis pipelie for BID-seq data" > /README.md

ENV PATH="/opt/pipeline/bin:/opt/miniconda/bin:$PATH"

# install system dependencies
RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y --no-install-recommends install tzdata apt-utils wget git make cmake xsltproc gcc g++ pkg-config zlib1g-dev python3 python3-distutils default-jre
# install pip
# wget https://bootstrap.pypa.io/get-pip.py -O /tmp/get-pip.py
# python3 /tmp/get-pip.py
# install miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.11.0-Linux-x86_64.sh -O /tmp/miniconda.sh && bash /tmp/miniconda.sh -b -p /opt/miniconda && rm -rf /tmp/miniconda.sh && export PATH="/opt/miniconda/bin:$PATH"
# install package by conda
RUN conda install mamba -n base -y -c conda-forge && mamba install -y -c conda-forge libzlib && mamba install -y -c bioconda falco bowtie2 star samtools==1.14 bedtools fastp seqtk blast subread numpy pandas pysam pyfaidx && mamba clean -afy && pip3 install snakemake==6.15.5 cutadapt==3.5 multiqc Flask
# install software directly
RUN wget https://github.com/getzlab/rnaseqc/releases/download/v2.4.2/rnaseqc.v2.4.2.linux.gz -O /opt/rnaseqc.gz && gunzip /opt/rnaseqc.gz && chmod +x /opt/rnaseqc
# install blast2bam
RUN git clone https://github.com/guyduche/Blast2Bam /opt/blast2bam && make -C /opt/blast2bam/ -j
# install rnaseqmut
RUN git clone https://github.com/davidliwei/rnaseqmut.git /opt/rnaseqmut && mkdir -p /opt/rnaseqmut/src/bamtools/build && cmake -B /opt/rnaseqmut/src/bamtools/build /opt/rnaseqmut/src/bamtools && make -C /opt/rnaseqmut/src/bamtools/build -j && make -C /opt/rnaseqmut/src -j
# install umicollapse
RUN git clone https://github.com/Daniel-Liu-c0deb0t/UMICollapse.git /opt/UMICollapse && mkdir -p /opt/UMICollapse/lib && wget -P /opt/UMICollapse/lib https://repo1.maven.org/maven2/com/github/samtools/htsjdk/2.19.0/htsjdk-2.19.0.jar && wget -P /opt/UMICollapse/lib https://repo1.maven.org/maven2/org/xerial/snappy/snappy-java/1.1.7.3/snappy-java-1.1.7.3.jar
# install cpup
RUN git clone https://github.com/y9c/cpup.git /opt/cpup && make -C /opt/cpup/ -j
RUN rm -rf /var/lib/apt/lists/* && apt-get purge -y wget git make cmake xsltproc gcc g++ pkg-config && apt-get clean

RUN mkdir -p /data
WORKDIR /data
ADD ./bin /opt/pipeline/bin
COPY ./Snakefile /opt/pipeline/Snakefile
COPY ./config.yaml /opt/pipeline/config.yaml
COPY ./bidseq /opt/pipeline/bin/bidseq
