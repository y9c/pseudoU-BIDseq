#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang <yech1990@gmail.com>
# Distributed under terms of the MIT license.
#
# Created: 2021-05-19 17:42

import gzip
import sys

import numpy as np

MAX_GAP_LEN = 5
UPSTREAM_NUM = 5
UPSTREAM_WEIGHT = 0.5


class Chunk:
    def __init__(self, samples):
        self.chrom = ""
        self.strand = ""
        self.counts = {sample: {} for sample in samples}
        self.nsample = len(samples)
        self.samples = samples
        self.pos = []


def process_count(count, strand, is_sort=False):
    """
    - count is modified
    """
    count_exteneded = []
    for pos, pos_stat in count.items():
        if pos_stat["D"]:
            for d_offset, d_length, d_count in pos_stat["D"]:

                def get_bases(pos_shift):
                    return [
                        count[d_pos]["base"]
                        for d_pos in range(
                            pos + 1 + pos_shift,
                            pos + d_length + pos_shift + 1,
                        )
                        if d_pos in count.keys()
                    ]

                bases = get_bases(0)
                pos_extended = [pos + d_offset]
                for pos_shift_sign in [-1, 1]:
                    pos_shift = 1
                    bases_shift = get_bases(pos_shift * pos_shift_sign)
                    while bases_shift and bases == bases_shift:
                        pos_extended.append(
                            pos + pos_shift * pos_shift_sign + d_offset
                        )
                        pos_shift += 1
                        bases_shift = get_bases(pos_shift * pos_shift_sign)
                count_exteneded.append(
                    [
                        tuple(np.arange(pos + 1, pos + 1 + d_length)),
                        tuple(sorted(pos_extended)),
                        d_count,
                    ]
                )
    if is_sort:
        count_exteneded = sorted(count_exteneded, key=lambda x: len(x[1]))
    return count_exteneded


def process_chunk(chunk):
    pos_max = max(chunk.pos)
    pos_min = min(chunk.pos)
    sort_counts = sorted(
        [
            [*c, sample]
            for sample, count in chunk.counts.items()
            for c in process_count(count, chunk.strand)
        ],
        key=lambda x: len(x[1]),
    )
    for gap_sites, gap_putatives, gap_count, sample in sort_counts:
        count = chunk.counts[sample]
        if len(gap_putatives) > 1:
            if chunk.strand == "+":
                gap_putatives = sorted(
                    sorted(gap_putatives, reverse=True),
                    key=lambda x: sum(
                        sum(
                            c[x + i]["gap"] * UPSTREAM_WEIGHT ** i
                            for c in chunk.counts.values()
                        )
                        for i in range(
                            min(pos_max - x + 1, UPSTREAM_NUM),
                        )
                    ),
                    reverse=True,
                )
            else:
                gap_putatives = sorted(
                    gap_putatives,
                    key=lambda x: sum(
                        sum(
                            c[x - i]["gap"] * UPSTREAM_WEIGHT ** i
                            for c in chunk.counts.values()
                        )
                        for i in range(
                            min(x - pos_min + 1, UPSTREAM_NUM),
                        )
                    ),
                    reverse=True,
                )
        gap_best = gap_putatives[0]
        if gap_best in count.keys():
            count[gap_best]["gap"] += gap_count
            for s in gap_sites:
                if s != gap_best:
                    if s in count.keys():
                        count[s]["depth"] += gap_count
                if gap_best not in gap_sites:
                    count[gap_best]["depth"] -= gap_count


def write_chunk(chunk):
    out = ""
    for p in chunk.pos:
        out += f"{chunk.chrom}\t{p}\t{chunk.strand}"
        for count in chunk.counts.values():
            out += f'\t{count[p]["depth"]}\t{count[p]["gap"]}'
        out += "\n"
    return out


if len(sys.argv) >= 3:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
else:
    input_file = "./test/test_input.tsv.gz"
    output_file = "./test/test_output.tsv.gz"


with gzip.open(input_file, "r") as f, gzip.open(output_file, "wb") as fo:
    names = [
        "depth",
        "A",
        "C",
        "G",
        "T",
        "N",
        "Gap",
        "Insert",
        "Delete",
        "istat",
        "dstat",
    ]
    _, _, _, _, *samples = f.readline().decode("utf-8").strip().split("\t")
    fo.write(
        (
            "\t".join(
                ["chr", "pos", "strand"]
                + [f"{s}_{t}" for s in samples for t in ["depth", "gap"]]
            )
            + "\n"
        ).encode()
    )
    pos_pre = 0
    chrom_pre = ""
    chunk = Chunk(samples)
    for l in f:
        chrom, pos, base, strand, *stats = (
            l.decode("utf-8").strip().split("\t")
        )
        pos = int(pos)
        if pos != pos_pre + 1 or chrom != chrom_pre:
            if len(chunk.pos) > 0:
                if any([c[p]["D"] for c in chunk.counts.values() for p in c]):
                    # do the analysis here
                    process_chunk(chunk)
                # write the output here
                fo.write(write_chunk(chunk).encode())
            chunk = Chunk(samples)
        chunk.pos.append(pos)
        for sample, stat in zip(samples, stats):
            *counts, istat, dstat = stat.split(",")
            pos_count = {x: int(y) for x, y in zip(names, counts)}
            pos_stat = {
                "base": base,
                "D": [],
                "depth": pos_count["A"]
                + pos_count["C"]
                + pos_count["G"]
                + pos_count["T"]
                + pos_count["N"],
                "gap": 0,
            }
            for d in dstat.split("|"):
                if d:
                    d_motif, d_count = d.split(":")
                    if len(d_motif) <= MAX_GAP_LEN:
                        if strand == "+":
                            d_offset = len(d_motif)
                        else:
                            d_offset = 1
                        pos_stat["D"].append(
                            (d_offset, len(d_motif), int(d_count))
                        )
            chunk.counts[sample][pos] = pos_stat
        chunk.chrom = chrom
        chunk.strand = strand
        chrom_pre = chrom
        pos_pre = pos

    if len(chunk.pos) > 0:
        if any([c[p]["D"] for c in chunk.counts.values() for p in c]):
            # do the analysis here
            process_chunk(chunk)
        # write the output here
        fo.write(write_chunk(chunk).encode())
