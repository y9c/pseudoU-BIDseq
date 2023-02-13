#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-02-10 22:41

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
        CALIBRATION_PARAMS[motif] = (float(c), float(s), float(b))


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
    if m not in CALIBRATION_PARAMS:
        return np.nan
    c, s, b = CALIBRATION_PARAMS[m]
    if x <= b:
        return 0.0
    elif x >= (c * (1 - s) + (1 - c) * b) / (1 - c * s):
        return 1.0
    else:
        return (b - x) / (c * (b + s - s * x - 1))


with open(input_file, "r") as fi, open(output_file, "w") as fo:
    header = fi.readline()
    fo.write(header)
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

    for line in fi:
        n_passed = 0
        records = line.strip().split("\t")
        motif = get_sequence(*records[:3])
        for _, idx_list in group_index.items():
            di, gi, dt, gt = [
                sum([int(records[i]) for i in idx]) for idx in idx_list
            ]
            if (
                dt >= GROUP_FILTER["min_treated_depth"]
                and di >= GROUP_FILTER["min_input_depth"]
                and gt >= GROUP_FILTER["min_treated_gap"]
                and gt / dt >= GROUP_FILTER["min_treated_ratio"]
                and (gt / dt) >= GROUP_FILTER["min_fold_ratio"] * (gi / di)
                and fit_calibration(gt / dt, motif)
                > GROUP_FILTER["min_treated_fraction"]
                and fisher_exact(
                    [[max(di - gi, 0), gi], [max(dt - gt, 0), gt]],
                    alternative="greater",
                ).pvalue
                < GROUP_FILTER["max_p_value"]
            ):
                n_passed += 1
            if n_passed >= GROUP_FILTER["min_passed_group"]:
                fo.write(line)
                break
