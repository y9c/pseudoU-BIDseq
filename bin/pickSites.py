#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-02-10 22:41

import gzip
import sys

import numpy as np
import pyfaidx
from scipy.stats import fisher_exact

GROUP_FILTER = snakemake.params.group_filter
GROUP_META = snakemake.params.group_meta
CALIBRATION_CURVE = snakemake.params.calibration_curve
REF_FASTA = snakemake.params.ref_fasta
input_file = snakemake.input[0]
output_file = snakemake.output[0]


CALIBRATION_PARAMS = {}
with open(CALIBRATION_CURVE, "r") as f:
    next(f)
    for line in f:
        motif, c, s, b = line.strip().split("\t")
        CALIBRATION_PARAMS[motif.replace("-", "T")] = (
            float(c),
            float(s),
            float(b),
        )


REF = pyfaidx.Fasta(REF_FASTA)


# assum the first 3 columns in the meta info of position
def get_sequence(chrom, pos, strand):
    pad = 2
    chrom = str(chrom)
    pos = int(pos)
    s5 = REF[chrom][pos - 1 - pad : pos - 1]
    s0 = REF[chrom][pos - 1]
    s3 = REF[chrom][pos : pos + pad]
    if strand == "+":
        return (
            s5.seq.rjust(pad, "N") + s0.seq + s3.seq.ljust(pad, "N")
        ).upper()
    elif strand == "-":
        return (
            s3.reverse.complement.seq.rjust(pad, "N")
            + s0.reverse.complement.seq
            + s5.reverse.complement.seq.ljust(pad, "N")
        ).upper()
    else:
        raise ValueError("Strand should be + or -")


def fit_calibration(x: float, m: str) -> float:
    #  m = m[len(m) // 2 - 2 : len(m) // 2 + 3]
    if np.isnan(x) or m not in CALIBRATION_PARAMS:
        return np.nan
    c, s, b = CALIBRATION_PARAMS[m]
    if x <= b:
        return 0.0
    elif x >= (c * (1 - s) + (1 - c) * b) / (1 - c * s):
        return 1.0
    else:
        return (b - x) / (c * (b + s - s * x - 1))


def open_file(filename, mode):
    if filename == "-":
        if "r" in mode:
            return sys.stdin
        elif "w" in mode:
            return sys.stdout
    elif filename.endswith(".gz"):
        if "r" in mode:
            return gzip.open(filename, "rt")
        elif "w" in mode:
            return gzip.open(filename, "wt")
    else:
        return open(filename, mode)


with open_file(input_file, "r") as fi, open_file(output_file, "w") as fo:
    header = fi.readline()
    header = header.strip().split("\t")
    # input depth, input gap, treated depth, treated gap
    group_index = {
        k: [
            [header.index(n + "_" + t) for n in v[l]]
            for l in ["input", "treated"]
            for t in ["depth", "gap"]
        ]
        for k, v in GROUP_META.items()
    }

    if GROUP_FILTER["combine_group_input"]:
        combined_input_depth = [
            header.index(n + "_depth")
            for v in GROUP_META.values()
            for n in v["input"]
        ]
        combined_input_gap = [
            header.index(n + "_gap")
            for v in GROUP_META.values()
            for n in v["input"]
        ]
        for k in GROUP_META:
            group_index[k][0] = combined_input_depth
            group_index[k][1] = combined_input_gap

    fo.write(
        "\t".join(
            header
            + [
                g + "_" + t
                for t in ["ratio", "fraction", "passed"]
                for g in group_index
            ]
        )
        + "\n"
    )

    for line in fi:
        n_passed = 0
        records = line.strip().split("\t")
        ratios = []
        fracs = []
        passeds = []
        motif = get_sequence(*records[:3])
        for g, idx_list in group_index.items():
            di, gi, dt, gt = [
                sum([int(records[i]) for i in idx]) for idx in idx_list
            ]
            ratio = gt / dt if dt > 0 else np.nan
            ratios.append(f"{ratio:.3f}" if not np.isnan(ratio) else "NaN")
            frac = fit_calibration(ratio, motif)
            fracs.append(f"{frac:.3f}" if not np.isnan(frac) else "NaN")
            if (
                dt >= GROUP_FILTER["min_treated_depth"]
                and di >= GROUP_FILTER["min_input_depth"]
                and gt >= GROUP_FILTER["min_treated_gap"]
                and ratio >= GROUP_FILTER["min_treated_ratio"]
                and ratio >= GROUP_FILTER["min_fold_ratio"] * (gi / di)
                and frac > GROUP_FILTER["min_treated_fraction"]
                and fisher_exact(
                    [[max(di - gi, 0), gi], [max(dt - gt, 0), gt]],
                    alternative="greater",
                ).pvalue
                < GROUP_FILTER["max_p_value"]
            ):
                passeds.append("1")
                n_passed += 1
            else:
                passeds.append("0")

        if n_passed >= GROUP_FILTER["min_passed_group"]:
            fo.write("\t".join(records + ratios + fracs + passeds) + "\n")
