#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2023 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2023-02-10 22:41


from scipy.stats import fisher_exact

GROUP_FILTER = snakemake.params.group_filter
GROUP_META = snakemake.params.group_meta
input_file = snakemake.input[0]
output_file = snakemake.output[0]


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
        records = line.strip().split("\t")
        passed = False
        n_passed = 0
        for g, (di_index, gi_index, dt_index, gt_index) in group_index.items():
            di = sum([int(records[i]) for i in di_index])
            gi = sum([int(records[i]) for i in gi_index])
            dt = sum([int(records[i]) for i in dt_index])
            gt = sum([int(records[i]) for i in gt_index])
            # min_treated_depth: 20
            # min_input_depth: 20
            # min_treated_gap: 5
            # min_treated_ratio: 0.02
            # min_fold_ratio: 2
            # max_p_value: 0.0001
            if (
                dt >= GROUP_FILTER["min_treated_depth"]
                and di >= GROUP_FILTER["min_input_depth"]
                and gt >= GROUP_FILTER["min_treated_gap"]
                and gt / dt >= GROUP_FILTER["min_treated_ratio"]
                and (gt / dt) >= GROUP_FILTER["min_fold_ratio"] * (gi / di)
                and fisher_exact(
                    [[max(di - gi, 0), gi], [max(dt - gt, 0), gt]],
                    alternative="greater",
                ).pvalue
                < GROUP_FILTER["max_p_value"]
            ):
                n_passed += 1
            if n_passed >= GROUP_FILTER["min_passed_group"]:
                fo.write(line)
                continue
