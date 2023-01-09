---
title: Algorithms
nav_exclude: false
nav_order: 5
---

<!-- prettier-ignore-start -->
# Algorithms
{: .fs-9 }
<!-- prettier-ignore-end -->

## realign gaps

Coming soon
{: .label .label-yellow }

`bowtie2` or other aligners use seed mapping to speed up, but also have some drawback on accurary.
Such as the example bellow, bowtie2 can not deal the gap (deletion) and the flanking sequence correctly.

- EXAMPLE 1: it dose not place the gap in the correct position and report a mismatch in the same time.

```
QRY: CGACG-TTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
     ||||. ||||||||||||||||||||||||||||||||||||
REF: CGACTGTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
```

After realign:

```
QRY: CGAC-GTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
     |||| |||||||||||||||||||||||||||||||||||||
REF: CGACTGTTTAATTAAAACAAAGCATCGCGAAGGCCCGCGGCG
```

- EXAMPLE 2: it treat gap region near the terminal as mutation.

```
QRY:    AGTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
        . |||||||||||||||||||||||||||||||||||||||||||||||
REF: (A)GTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
```

After realign:

```
QRY: AG-TGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
     || |||||||||||||||||||||||||||||||||||||||||||||||
REF: AGTTGAAAAGAACTTTGAAGAGAGAGTTCAAGAGGGCGTGAAACCGTTAA
```

## Adjust gaps

Coming soon
{: .label .label-yellow }
