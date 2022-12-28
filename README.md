[![Docker](https://img.shields.io/docker/pulls/y9ch/bidseq.svg)](https://hub.docker.com/r/y9ch/bidseq)

# &Psi;-BID-seq

## Overview of the workflow

<p align="center">
  <a href="https://y9c.github.io/pseudoU-BIDseq/Overall-Workflow#gh-light-mode-only">
    <img src="./docs/scheme.svg" />
  </a>
  <a href="https://y9c.github.io/pseudoU-BIDseq/Overall-Workflow#gh-dark-mode-only">
    <img src="./docs/scheme_dark.svg" />
  </a>
</p>

## How to use?

A [docker image](https://hub.docker.com/r/y9ch/bidseq) containing the source code and dependencies has been published for reproducibility. You can run it using the [singularity](https://sylabs.io/singularity) container runtime.

The entire analysis can be completed in just three steps:

1. **Specific the path (with label) of both rawdata and references for your project in a YAML format.**

   <details>
     <summary><code>data.yaml</code> for example<sup>(Click to expand)</sup></summary>

   ```yaml
   reference:
     contamination:
       fa: ../ref/contamination.fa
       bt2: ../ref/contamination
     genes:
       fa: ../ref/genes.fa
       fai: ../ref/genes.fa.fai
       bt2: ../ref/genes
     genome:
       fa: /data/reference/genome/Mus_musculus/GRCm39.fa
       fai: /data/reference/genome/Mus_musculus/GRCm39.fa.fai
       star: /data/reference/genome/Mus_musculus/star/GRCm39.release108

   samples:
     mESCWT-rep1-input:
       data:
         - R1: ../test/IP16.fastq.gz
       group: mESCWT
       treated: false
     mESCWT-rep2-input:
       data:
         - R1: ../test/IP17.fastq.gz
       group: mESCWT
       treated: false
     mESCWT-rep1-treated:
       data:
         - R1: ../test/IP4.fastq.gz
       group: mESCWT
       treated: true
     mESCWT-rep2-treated:
       data:
         - R1: ../test/IP5.fastq.gz
       group: mESCWT
       treated: true
   ```

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Run-the-pipeline.html#refer-rawdata-and-references-in-the-configuration-file) on how to customize._

   </details>

2. **Run all the analysis by one command**:

   ```bash
   singularity exec docker://y9ch/bidseq:latest bidseq
   ```

    <details>
      <summary>default settings<sup>(Click to expand)</sup></summary>

   - default config file: `data.yaml`
   - default output dir: `./results`
   - default jobs in parallel: `48`

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Run-the-pipeline.html#customized-analysis-parameters) on how to customize._

   </details>

3. **View the analytics report and use the &Psi; sites for downstream analysis**.

   The output of all the steps will be in the specific working directory.

## Documentation

...

(**⚠ Not yet published.**)

## Citation

- Qing Dai†_, Li-Sheng Zhang†_, Hui-Lung Sun†, Kinga Pajdzik†, Lei Yang, Chang Ye, Cheng-Wei Ju, Shun Liu, Yuru Wang, Zhong Zheng, Linda Zhang, Bryan Harada, Xiaoyang Dou, Iryna Irkliyenko, Xinran Feng, Wen Zhang, Tao Pan, Chuan He\*. Quantitative sequencing using BID-seq uncovers abundant pseudouridines in mammalian mRNA at base resolution. Nat Biotechnol (2022). https://doi.org/10.1038/s41587-022-01505-w

&nbsp;

<p align="center">
  <img
    src="https://raw.githubusercontent.com/catppuccin/catppuccin/dev/assets/footers/gray0_ctp_on_line.svg?sanitize=true"
  />
</p>
<p align="center">
  Copyright &copy; 2021-present
  <a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
  <a href="https://github.com/y9c/pseudoU-BIDseq/blob/master/LICENSE"
    ><img
      src="https://img.shields.io/static/v1.svg?style=for-the-badge&label=License&message=GPLv3&logoColor=d9e0ee&colorA=282a36&colorB=c678dd"
  /></a>
</p>
