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

   You can also copy and edit from this [template](https://github.com/y9c/pseudoU-BIDseq/blob/main/test/data.yaml).

   _Read the [more details](https://y9c.github.io/pseudoU-BIDseq/Step-by-step-instruction.html#define-settings-in-the-configure-file) on how to customize._

2. **Run all the analysis by one command**:

   ```bash
   apptainer run docker://y9ch/bidseq
   ```

   default
   {: .label .label-blue }
   The pipeline will load configure file named `data.yaml` under the current directory.

   {: .note }

   > How to run apptainer on computation nodes without internet acess?
   >
   > 1. (On the login node with internet connection) Run `module load apptainer` to mount the apptainer utils, if it is not installed by default.
   >
   > 2. (On the login node with internet connection) Build the `bidseq_latest.sif` file using the command `apptainer pull docker://y9ch/bidseq`.
   >
   > 3. (On the computation node) Run `apptainer run bidseq_lastest.sif -c data.yaml` to start the pipeline.
   >    Note that most HPC systems mount directories in a complex manner. Therefore, you need to find out the actual path by executing `realpath ./` and specify this output into `apptainer` using `apptainer run -B /the/real/path ...`

   {: .note }

   > If your configure file is not named as `data.yaml`, add `-c your_file_name.yaml` arg after the command to customize.

3. **View the analytics reports and filtered sites.**

   default
   {: .label .label-blue }
   3 folder will be created in the working directory (default: `workspace`),

   {: .note }

   > - trimming, mapping, deduping reports are in `report_reads` folder, with key numbers in all the steps reported in one webpage<sup>([example](https://y9c.github.io/pseudoU-BIDseq/readsStats))</sup>.
   > - filtered sites for &Psi; sites detection are in `filter_sites` folder. These sites are only passed the _simplest filtering_, you can apply customized threshold into them based your data type and quality.
   > - processed mapping results (_.bam_) are in `align_bam` folder. You can zoom into location that you interested in IGV.
