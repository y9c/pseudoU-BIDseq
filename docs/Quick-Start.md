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

   For example, write down and save the following block into a text file and named it as `data.yaml`.

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

   You can also copy and edit from this [template](test/data.yaml).

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Step-by-step-instruction.html#define-settings-in-configure-file) on how to customize._

2. **Run all the analysis by one command**:

   ```bash
   apptainer run docker://y9ch/bidseq:v1
   ```

   default
   {: .label .label-blue }
   The pipeline will load configure file named as `data.yaml` under current directory

   {: .note }

   > - Customized configure file with `-c` argument. (default: `data.yaml`)
   > - Customized number of jobs/cores in parallel `-j` argument. (default: `48`)

   _Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Step-by-step-instruction.html#customized-settings-in-command-line) on how to customize._

3. **View the analytics report and pileup table.**

   default
   {: .label .label-blue }
   3 folder are will be created in the working directory (default: `workspace`),

   {: .note }

   > ├── align_bam
   > ├── **report_reads**
   > └── **filter_sites**
   >
   > - trimming, mapping, deduping reports are in `report_reads` folder, with key numbers in all the steps reported in one webpage<sup>([example](https://y9c.github.io/pseudoU-BIDseq/readsStats))</sup>.
   > - deleted sites for &Psi; sites detection are in `filter_sites` folder. These sites are only passed the _simplest filtering_, you can apply customized threadfolds into them based your data type and quality.
   > - processed mapping results (_.bam_) are in `align_bam` folder. You can zoom into location that you interested in IGV.
