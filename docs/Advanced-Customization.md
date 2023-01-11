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

According to the fragment size of BID-seq libraries, SE100 or SE150 sequencing mode is sufficient for most of the RNA fragments.
To reduce cost, you do not need to run pair-end sequencing for most of the samples.
But if you have a longer insert fragment size, or want to improve sequencing quality by PE mode, this pipeline also supports analyzing data from PE mode.

And the setup is as simple as adding a single line into the YAML configure file.

```yaml
samples:
  mESCWT-rep1-input:
    data:
      - R1: ./test/IP16_R1.fastq.gz
        R2: ./test/IP16_R2.fastq.gz
```

`R2: ./test/IP16_R2.fastq.gz` is added under the data tag.

{: .note }
Read 2 file is labeled with an `R2` tag, and it is in the same indent as the `R1` tag. There is no hyphen symbol before the `R2` tag, because R2 and R1 are in pairs.

## Multiple sequencing runs

If you sequenced the same library in multiple sequencing flowcells or added more reads for your sample, you do not need to combine the data before running this pipeline.
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
Sequencing data is combined after read alignment, rather than at the first step of the analysis. This strategy can save computation resources and energy. For example, sometimes you run the sequencing for your libraries but found that the data is not sufficient after the analysis.
You then add extra sequencing data for this library. In this pipeline, only newly generated data need to be aligned.

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

## Customized adapter / inline barcode

You can customize the adapter sequencing if you are not using the adapter (\_NNNNN_AGATCGGAAGAGCACACGTCT) provided by the protocol.

By default, only 5 N are added on the 3' adapter, which is used as an inline barcode. But it is also possible to add an inline barcode as the one used in _Nature__ Biotech.__ [paper](https://www.nature.com/articles/s41587-022-01505-w#Sec12) (_NNNNN_<u>ATCACG</u>AGATCGGAAGAGCACACGTCT).

There are two ways to specify the inline barcode.

- You can apply the global setting for all the samples in the configure file by adding:

```yaml
barcode: ATCACG
```

- If only some of the libraries are with inline barcodes, while others are not, you can specific inline barcodes for each sample respectively. Leave it blank (default) means without an inline barcode.

```yaml
samples:
  mESCWT-rep1-input:
    data:
      - R1: ./test/IP16_sequencing_run1.fastq.gz
    barcode: ATCACG
```

## Customized cutoff for pre-filtering

Not recommended
{: .label .label-red }

Add the following **whole** block into the configure file (`data.yaml` for example), and adjust the parameters.

```yaml
cutoff:
  min_match_prop: 0.8
  min_group_gap: 5
  min_group_depth: 10
  min_group_ratio: 0.01
  min_group_num: 1
```

- `min_match_prop`: STAR `outFilterMatchNminOverLread`
- `min_group_gap`: prefilter sites show >= x gaps in total among all samples in each group
- `min_group_depth`: prefilter sites show >= x sequencing coverage in total among all samples in each group
- `min_group_ratio`: prefilter sites show >= x deletion ratio among all samples in each group
- `min_group_num`: only analysis putative sites that show pass prefilter in >= x groups

## Cache internal files to speed up

Add the following setting in the configure file to turn on `keep_internal` (default: false). Once internal files, including reference index and mapping bam files, are cached, you do not need to re-run some steps of the pipeline when you add more sequencing reads.

## Keep discarded reads for debugging purposes.

```yaml
keep_discarded: true
```

Once this parameter is set as true (default: False), untrimmed, too-short, unmapped... reads will be saved.
_Most of the time, you do not need these files. Set it as `false` to save storage._
