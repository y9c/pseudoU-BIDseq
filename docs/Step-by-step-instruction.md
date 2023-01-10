---
title: Step by Step Instruction
nav_exclude: false
nav_order: 3
math: mathjax3
---

<!-- prettier-ignore-start -->
# Step by Step Instruction
{: .fs-9 }
<!-- prettier-ignore-end -->

## Define settings in configure file

> Important parameters of the `data.yaml` file.

| --------- | ------------- | ------- | --------------------    | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| --------- | ------------- | ------- | ----------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Parameter |               |         | Default value (type)    | Description                                                                                                                                                                                                              |
| Level 1   | Level 2       | Level 3 | ^^                      | ^^                                                                                                                                                                                                                       |
| workdir   | -             | -       | workspace (string)      | (_optional_) The path to the working directory.                                                                                                                                                                          |
| tempdir   | -             | -       | workspace/.tmp (string) | (_optional_) The path to the temp file directory.                                                                                                                                                                        |
| cores     | -             | -       | 48 (int)                | (_optional_) Max number of cores for all the running jobs at the same time.                                                                                                                                              |
| reference | contamination | fa      | - (string)              | (_optional_) Reference file for putative contamination in the samples                                                                                                                                                    |
| ^^        | ^^            | bt2     | Auto (string)           | (_optional_) `bowtie2` index for contamination reference                                                                                                                                                                 |
| ^^        | genes         | fa      | - (string)              | (**required**) reference file for selected rRNA, tRNA and snoRNA genes                                                                                                                                                   |
| ^^        | ^^            | bt2     | Auto (string)           | (_optional_) `bowtie2` index for selected genes                                                                                                                                                                          |
| ^^        | genome        | fa      | - (string)              | (**required**) reference genome for the species you study                                                                                                                                                                |
| ^^        | ^^            | star    | - (string)              | (**required**) `STAR` index for reference genome                                                                                                                                                                         |
| sample    | (sample name) | data    | - (list)                | (**required**) list of files path for sequencing reads in fastq format. Multiple sequencing run for one sample is supported. You do not need to combine the input fastq files before passing the data into the pipeline. |
| ^^        | ^^            | group   | - (string)              | (**required**) group ID for this sample                                                                                                                                                                                  |
| ^^        | ^^            | treated | true (boolean)          | (_optional_) true for BS treated sample, false for untreated control sample.                                                                                                                                             |
| ^^        | (…)           | (...)   | (…)                     | (_optional_) other samples. Adding any number of entries is supported.                                                                                                                                                   |

- `workdir`
- `tempdir`
- `cores`
- ...

_Read the [documentation](https://y9c.github.io/pseudoU-BIDseq/Advanced-Customization.html) on how to customize._

## Customized settings in command line

- Customized configure file with `-c` argument. (default: `data.yaml`)
- Customized number of jobs/cores in parallel `-j` argument. (default: `48`)

## Post-filter &Psi; sites

### The reliablity ($p$) of the &Psi;-modified sites

A p-value from **One-sided Fisher’s Exact Test** can be used for evaluating the reliabity of each site.

As the table bellow,

- $\sum{d_{i}}$ is total number of coverage in input libraries
- $\sum{d_{t}}$ is total number of coverage in treated libraries
- $\sum{g_{i}}$ is total number of gaps in input libraries
- $\sum{g_{t}}$ is total number of gaps in treated libraries

| --------- | -----------------------     | -----------------           |
| --------- | --------------------------- | --------------------------- |
|           | Input                       | Treated                     |
| gap       | $\sum{g_{i}}$               | $\sum{g_{t}}$               |
| T         | $\sum{d_{i}} - \sum{g_{i}}$ | $\sum{d_{t}} - \sum{g_{t}}$ |

$$
\displaystyle p={\frac{  \displaystyle{ {\sum{d_{i}}}\choose{\sum{g_{i}}} } \displaystyle{ {\sum{d_{t}}}\choose{\sum{g_{t}}} }  }{\displaystyle{ {\sum{d_{i}}+\sum{d_{t}}}\choose{\sum{g_{i}}+\sum{g_{t}}} }}}
$$

{: .note }

> You can add p.value to the table generated by this pipeline with a R script.
>
> Some example rows from the file:
>
> ```
> chr	pos	strand	mESCWT-rep1-input_depth	mESCWT-rep1-input_gap	mESCWT-rep2-input_depth	mESCWT-rep2-input_gap	mESCWT-rep3-input_depth	mESCWT-rep3-input_gap	mESCWT-rep1-treated_depth	mESCWT-rep1-treated_gap	mESCWT-rep2-treated_depth	mESCWT-rep2-treated_gap	mESCWT-rep3-treated_depth	mESCWT-rep3-treated_gap
> rRNA_Rn45s	43	+	112	28	74	25	72	25	108	25	96	26	88	22
> rRNA_Rn45s	4913	+	416	1	362	0	396	0	203	29	216	29	196	26
> rRNA_Rn45s	10032	+	182	0	150	2	161	0	179	0	173	0	157	0
> rRNA_Rn45s	10970	+	124	2	118	2	134	5	142	5	133	2	113	4
> rRNA_Rn45s	11491	+	327	0	315	2	278	0	153	129	159	130	158	128
> ...
> ```
>
> Example code snippet:
>
> ```R
> # customized your file name and sample name in these 4 rows
> fn_in <- "sites.tsv.gz" # input file anme
> fn_out <- "sites_pvalue.tsv.gz" # output file name
> col_i <- c("mESCWT-rep1-input", "mESCWT-rep2-input", "mESCWT-rep3-input") # sample of input libs
> col_t <- c("mESCWT-rep1-treated", "mESCWT-rep2-treated", "mESCWT-rep3-treated") # sample of treated libs
>
> # calculate p.value
> get_p <- function(r) {
>   fisher.test(as.table(matrix(as.numeric(r[c(1:4)]), ncol = 2)), alt = "greater")$p.value
> }
> df <- read.table(gzfile(fn_in), sep = "\t", header = T, check.names = F)
> di <- rowSums(df[paste0(col_i, "_depth")])
> gi <- rowSums(df[paste0(col_i, "_gap")])
> dt <- rowSums(df[paste0(col_t, "_depth")])
> gt <- rowSums(df[paste0(col_t, "_gap")])
> df$p.value <- apply(cbind(di - gi, gi, dt - gt, gt), 1, get_p)
> write.table(df, gzfile(fn_out), sep = "\t", row.names = F, quote = F)
> ```

### The stoichiometry ($f$) of the &Psi;-modified sites

$f$ can be precisely calculated by applying the deletion ratio ($r = \displaystyle{\frac{g}{d}}$, where $g$ is the deletion number and $d$ is the sequencing coverage) to the calibration curves.

$$
f=\frac{b - r}{c \times (b + s - s \times r - 1)}
$$

, where $c$ is the <u>conversion ratio</u>, $s$ is the <u>RT dropout proportion</u> and $b$ is the <u>background noise</u>.
