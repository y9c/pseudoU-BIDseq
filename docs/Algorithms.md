---
title: Algorithms
nav_exclude: false
nav_order: 5
---

<!-- prettier-ignore-start -->
# Algorithms
Overall Workflow
{: .fs-9 }
<!-- prettier-ignore-end -->

## realign gaps

bowtie2 aligner which is based on Burrows-Wheeler Transform algorithm would report the "best" alignmet each time.
Such as the example bellow. bowtie2 can not align the gap (deletion) correctly, and report a mismatch next to the gap postion in the same time.

```
A00639:852:HH5Y5DRXY:1:2103:10836:12179
CGACG-TTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
||||. ||||||||||||||||||||||||||||||||||||
CGACTGTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
```

This can be fixed by gap realignment, and have a better result:

```
A00639:852:HH5Y5DRXY:1:2103:10836:12179
CGAC-GTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
|||| |||||||||||||||||||||||||||||||||||||
CGACTGTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
```

mismatch at the terminal of sequence.

```
A00639:852:HH5Y5DRXY:1:2101:24876:34554
   AGTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
   . |||||||||||||||||||||||||||||||||||||||||||||||
(A)GTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
```

into:

```
A00639:852:HH5Y5DRXY:1:2101:24876:34554
AG-TGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
|| |||||||||||||||||||||||||||||||||||||||||||||||
AGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
```
