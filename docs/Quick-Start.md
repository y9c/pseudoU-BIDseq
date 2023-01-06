---
title: Quick Start
nav_exclude: false
nav_order: 2
---

<!-- prettier-ignore-start -->
# Quick Start
{: .fs-9 }
<!-- prettier-ignore-end -->

1. **Specific the path of references (_.fasta_) and samples (_.fastq_) in a configure file (_.YAML_).**

   <details>
     <summary><code>data.yaml</code> for example<sup>(Click to expand)</sup></summary>

   ```yaml
   reference:
     contamination:
       fa: ./ref/contamination.fa
     genes:
       fa: ./ref/genes.fa
     genome:
       fa: /data/reference/genome/Mus_musculus/GRCm39.fa
       star: /data/reference/genome/Mus_musculus/star/GRCm39.release108

   samples:
     mESCWT-rep1-input:
       data:
         - R1: ./test/IP16.fastq.gz
       group: mESCWT
       treated: false
     mESCWT-rep1-treated:
       data:
         - R1: ./test/IP4.fastq.gz
       group: mESCWT
       treated: true
     mESCWT-rep2-treated:
       data:
         - R1: ./test/IP5.fastq.gz
       group: mESCWT
       treated: true
   ```

   You can copy and edit from this [template](test/data.yaml).

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Run-the-pipeline.html#refer-rawdata-and-references-in-the-configuration-file) on how to customize._

   </details>

2. **Run all the analysis by one command**:

   ```bash
   apptainer run docker://y9ch/bidseq:v1
   ```

    <details>
      <summary>default settings<sup>(Click to expand)</sup></summary>

   default
   {: .label .label-blue }

   - default config file: `data.yaml`
   - default output dir: `./workspace`
   - default jobs in parallel: `48`

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Run-the-pipeline.html#customized-analysis-parameters) on how to customize._

   </details>

3. **View the analytics report and pileup table.**

    <details>
      <summary>3 folder are will be created in the working directory.<sup>(Click to expand)</sup></summary>

   default
   {: .label .label-blue }

   - trimming, mapping, deduping reports are in `report_reads` folder, with key numbers in all the steps reported in one webpage<sup>([example](https://y9c.github.io/pseudoU-BIDseq/readsStats))</sup>.
   - deleted sites for &Psi; sites detection are in `pileup_filtered` folder. These sites are only passed the _simplest filtering_, you can apply customized threadfolds into them based your data type and quality.
   - processed mapping results (_.bam_) are in `drop_duplicates` folder. You can zoom into location that you interested in IGV.
   </details>
