---
title: Advanced Customization
nav_exclude: false
nav_order: 4
---

<!-- prettier-ignore-start -->
# Advanced Customization
{: .fs-9 }
<!-- prettier-ignore-end -->

## Pair-end mode

According to the size of BID-seq libraries, SE100 or SE150 sequencing mode is sufficent for most of the RNA fragments.
To reduce cost, you do not need to run pair-end sequencing for most of the samples.
But if you have longer insert fragment size, or want to imporve better sequencing quality by PE mode, this pipeline also support analyzing data from PE mode.

And the set up is as simple as adding a single line to into the YAML configure file.

```yaml
samples:
  mESCWT-rep1-input:
    data:
      - R1: ./test/IP16_R1.fastq.gz
        R2: ./test/IP16_R2.fastq.gz
```

`R2: ./test/IP16_R2.fastq.gz` is added under the data tag.

{: .note }
Read 2 file is labeled with a `R2:` tag, and it is in the same indent as the `R1` tag. There is no hypen symbol before `R2` tag, because R2 and R1 are in pairs.

## multiple sequencing runs

If you sequenced the same library in multiple sequencing flowcells, or added more reads for your sample, you do not need to combined the data before running this pipeline.
You can add multiple runs under the data tag for one sample. When you start the pipeline, it will combine the data for this sample automatically.

```yaml
samples:
  mESCWT-rep1-input:
    data:
      - R1: ./test/IP16_sequencing_run1.fastq.gz
      - R1: ./test/IP16_sequencing_run2.fastq.gz
      - R1: ./test/IP16_sequencing_run2.fastq.gz
```

{: .note }
Sequencing data is combined after read alignment, rather than at the first step of the analysis. This stategy can save computation resource and energy. For example, sometime you run the sequencing for your libraries, but found that the data is not sufficent after the analyze.
You then add extra sequencing data for this library. In this pipeline, only new generated data need to be aligned.

## Use pre-analyzed bam file for &Psi; sites detection only

```yaml
samples:
  mESCWT-rep1-treated:
    bam:
      gene: A1.gene.bam
      genome: A1.genome.bam
  mESCKO-rep1-treated:
    bam:
      gene: B1.gene.bam
      genome: B1.genome.bam
```

## Customized adapter (inline barcode)

If you customized your adapter sequencing if you are not using the adpter provided in the protocol.

Coming soon
{: .label .label-yellow }
Customized

## Customized path of tools

Coming soon
{: .label .label-yellow }
Customized

## Cache internal files to speed up

Add the following setting in the configure file, to turn on `keep_internal` (default: false). Once internal files, including refernce index and mapping bam files, are cached, you do not need to re-run some analysis when you add more sequening reads.

```yaml
keep_internal: true
```

## keep dicarded reads for debugging purpose.

```yaml
keep_discarded: true
```

Once this parameter is true (default: False), untrimmed, too-short, unmapped... reads will be saved.
_Most of the time, you do not need these files. Set it as `false` to save storage._

## Mics

Coming soon
{: .label .label-yellow }

```yaml
trim_p5: false
greedy_mapping: false
```
