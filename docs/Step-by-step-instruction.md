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

## define settings in configure file

> Important parameters of the `data.yaml` file.

| --------- | ------------- | ------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| --------- | ------------- | ------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Parameter |               |         | Default value (type) | Description                                                                                                                                                                                                              |
| Level 1   | Level 2       | Level 3 | ^^                   | ^^                                                                                                                                                                                                                       |
| workdir   | -             | -       | workspace (string)   | (_optional_) The path to the working directory.                                                                                                                                                                          |
| cores     | -             | -       | 48 (int)             | (_optional_) Max number of cores for all the running jobs at the same time.                                                                                                                                              |
| reference | contamination | fa      | - (string)           | (_optional_) Reference file for putative contamination in the samples                                                                                                                                                    |
| ^^        | ^^            | bt2     | Auto (string)        | (_optional_) `bowtie2` index for contamination reference                                                                                                                                                                 |
| ^^        | genes         | fa      | - (string)           | (**required**) reference file for selected rRNA, tRNA and snoRNA genes                                                                                                                                                   |
| ^^        | ^^            | bt2     | Auto (string)        | (_optional_) `bowtie2` index for selected genes                                                                                                                                                                          |
| ^^        | genome        | fa      | - (string)           | (**required**) reference genome for the species you study                                                                                                                                                                |
| ^^        | ^^            | star    | - (string)           | (**required**) `STAR` index for reference genome                                                                                                                                                                         |
| sample    | (sample name) | data    | - (list)             | (**required**) list of files path for sequencing reads in fastq format. Multiple sequencing run for one sample is supported. You do not need to combine the input fastq files before passing the data into the pipeline. |
| ^^        | ^^            | group   | - (string)           | (**required**) group ID for this sample                                                                                                                                                                                  |
| ^^        | ^^            | treated | true (boolean)       | (_optional_) true for BS treated sample, false for untreated control sample.                                                                                                                                             |
| ^^        | (…)           | (...)   | (…)                  | (_optional_) other samples. Adding any number of entries is supported.                                                                                                                                                   |

## customized settings in command line

- Customized configure file with `-c` argument. (default: `data.yaml`)
- Customized number of jobs/cores in parallel `-j` argument. (default: `48`)

## post filter &Psi; sites

### The stoichiometry ($f$) of the &Psi;-modified sites

$f$ can be precisely calculated by applying the deletion ratio ($r = \displaystyle{\frac{g}{d}}$, where $g$ is deletion number and $d$ is sequencing coverage) to the calibration curves.

$$
f=\frac{b - r}{c \times (b + s - s \times r - 1)}
$$

, where $c$ is conversion ratio conversion ratio, $s$ is the RT dropout proportion and $b$ is the background noise.

### The reliablity ($p$) of the &Psi;-modified sites

The p-value from One-sided Fisher’s Exact Test for number of T bases in ‘BID-seq’ libraries ($\sum{d_{t}} - \sum{g_{t}}$), number of deletions in ‘BID-seq’ libraries ($\sum{g_{t}}$), number of T bases in ‘Input’ libraries ($\sum{d_{i}} - \sum{g_{i}}$), number of deletions in ‘Input’ libraries ($\sum{g_{i}}$) is less than $10^{-4}$.

| --------- | -----------------------     | -----------------           |
| --------- | --------------------------- | --------------------------- |
|           | Input                       | Treated                     |
| T         | $\sum{d_{i}} - \sum{g_{i}}$ | $\sum{d_{t}} - \sum{g_{t}}$ |
| gap       | $\sum{g_{i}}$               | $\sum{g_{t}}$               |

$$
\displaystyle p={\frac{  \displaystyle{ {\sum{d_{i}}}\choose{\sum{g_{i}}} } \displaystyle{ {\sum{d_{t}}}\choose{\sum{g_{t}}} }  }{\displaystyle{ {\sum{d_{i}}+\sum{d_{t}}}\choose{\sum{g_{i}}+\sum{g_{t}}} }}}
$$
