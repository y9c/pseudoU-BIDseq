[![Docker](https://img.shields.io/docker/pulls/y9ch/bidseq.svg)](https://hub.docker.com/r/y9ch/bidseq)

# &Psi;-SAC-seq

## Overview of the workflow

<p align="center">
  <a href="https://y9c.github.io/m6A-SACseq/Overall-Workflow#gh-light-mode-only">
    <img src="./docs/scheme.svg" />
  </a>
  <a href="https://y9c.github.io/m6A-SACseq/Overall-Workflow#gh-dark-mode-only">
    <img src="./docs/scheme_dark.svg" />
  </a>
</p>

## How to use?

A [docker image](https://hub.docker.com/r/y9ch/sacseq) containing the source code and dependencies has been published for reproducibility. You can run it using the [singularity](https://sylabs.io/singularity) container runtime.

The entire analysis can be completed in just three steps:

1. **Specific the path (with label) of both rawdata and references for your project in a YAML format.**

<details>
  <summary><code>data.yaml</code> for example<sup>(Click to expand)</sup></summary>

```yaml
samples:
  HeLa-WT:
    input:
      rep1:
        - R1: ./rawdata/HeLa-WT-polyA-input-rep1-run1_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-input-rep1-run1_R2.fq.gz
        - R1: ./rawdata/HeLa-WT-polyA-input-rep1-run2_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-input-rep1-run2_R2.fq.gz
      rep2:
        - R1: ./rawdata/HeLa-WT-polyA-input-rep2-run1_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-input-rep2-run1_R2.fq.gz
        - R1: ./rawdata/HeLa-WT-polyA-input-rep2-run2_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-input-rep2-run2_R2.fq.gz
    treated:
      rep1:
        - R1: ./rawdata/HeLa-WT-polyA-treated-rep1-run1_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-treated-rep1-run1_R2.fq.gz
        - R1: ./rawdata/HeLa-WT-polyA-treated-rep1-run2_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-treated-rep1-run2_R2.fq.gz
      rep2:
        - R1: ./rawdata/HeLa-WT-polyA-treated-rep2-run1_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-treated-rep2-run1_R2.fq.gz
        - R1: ./rawdata/HeLa-WT-polyA-treated-rep2-run2_R1.fq.gz
          R2: ./rawdata/HeLa-WT-polyA-treated-rep2-run2_R2.fq.gz
references:
  spike:
    fa: ./ref/spike-in.fa
    bt2: ./ref/spike-in.fa
  spikeN:
    blast: ./ref/spike-in_with_N
  rRNA:
    fa: ./ref/Homo_sapiens.GRCh38.rRNA.fa
    bt2: ./ref/Homo_sapiens.GRCh38.rRNA
  smallRNA:
    fa: ./ref/Homo_sapiens.GRCh38.smallRNA.fa
    bt2: ./ref/Homo_sapiens.GRCh38.smallRNA
  genome:
    fa: ./ref/Homo_sapiens.GRCh38.genome.fa
    star: ./ref/Homo_sapiens.GRCh38.genome
    gtf: ./ref/Homo_sapiens.GRCh38.genome.gtf
    fai: ./ref/Homo_sapiens.GRCh38.genome.fa.fai
    gtf_collapse: ./ref/Homo_sapiens.GRCh38.genome.collapse.gtf
  contamination:
    fa: ./ref/contamination.fa
    bt2: ./ref/contamination
```

_Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Step-by-step-instruction.html#refer-rawdata-and-references-in-the-configuration-file) on how to customize._

</details>

2. **Run all the analysis by one command**:

```bash
singularity exec docker://y9ch/sacseq:latest sacseq
```

<details>
  <summary>default settings<sup>(Click to expand)</sup></summary>

- default config file: `data.yaml`
- default output dir: `./results`
- default jobs in parallel: `48`

_Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Step-by-step-instruction.html#customized-analysis-parameters) on how to customize._

</details>

3. **View the analytics report and use the &Psi; sites for downstream analysis**.

The output of all the steps will be in one folder (`./results`) under the current path. A webpage report of all the analysis will be in `./results/report.html` <sup>([example](https://y9c.github.io/pseudoU-BIDseq/demo_output.html))</sup>.

## Documentation

https://y9c.github.io/m6A-SACseq/

## Citation

TBD
