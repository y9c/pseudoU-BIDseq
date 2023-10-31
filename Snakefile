import sys
from snakemake.utils import min_version
from collections import defaultdict

if sys.version_info < (3, 6):
    sys.exit("Python 3.6 or later is required.\n")
min_version("7.0")


WORKDIR = os.path.relpath(
    config.get("workdir", "workspace"), os.path.dirname(workflow.configfiles[-1])
)
TEMPDIR = os.path.relpath(config.get("tempdir", os.path.join(WORKDIR, ".tmp")), WORKDIR)
INTERNALDIR = "internal_files"


workdir: WORKDIR


CALI = config.get("calibration_curves", "")
CALI = (
    CALI
    if os.path.isabs(CALI)
    else os.path.expanduser(CALI)
    if CALI.startswith("~")
    else os.path.relpath(CALI, WORKDIR)
)

REF = config["reference"]
for k, v in REF.items():
    for k2, v2 in v.items():
        v2 = os.path.expanduser(v2)
        REF[k][k2] = v2 if os.path.isabs(v2) else os.path.relpath(v2, WORKDIR)


def parse_barcode(b):
    # NNNNNXXX-XXXNNNNNATCACG
    if "-" in b:
        b1m1, m2b2i = b.split("-")
    else:
        b1m1, m2b2i = "", b
    # parse left
    m1 = b1m1.lstrip("N")
    b1 = b1m1[: len(b1m1) - len(m1)]
    if not all([x == "X" for x in m1]):
        raise ValueError(f"5'-mask {m1} is not all N.")
    # parse right
    b2i = m2b2i.lstrip("X")
    m2 = m2b2i[: len(m2b2i) - len(b2i)]
    i = b2i.lstrip("N")
    b2 = b2i[: len(b2i) - len(i)]
    # check inline barcode is "ATGC"
    for x in i:
        if x not in "ATGCatgc":
            raise ValueError(f"Inline barcode {i} is not all A/T/G/C.")
    return {
        "inline": i,
        "umi5": len(b1),
        "umi3": len(b2),
        "mask5": len(m1),
        "mask3": len(m2),
    }


REFTYPE = ["genes", "genome"]
GROUP2SAMPLE = defaultdict(lambda: defaultdict(list))
SAMPLE_IDS = []
SAMPLE2RUN = defaultdict(dict)
SAMPLE2BARCODE = defaultdict(dict)
SAMPLE2BAM = defaultdict(dict)
for s, v2 in config["samples"].items():
    SAMPLE_IDS.append(s)
    if v2.get("treated", True):
        if isinstance(v2["group"], list):
            raise ValueError("treated samples can only be in one group")
        GROUP2SAMPLE[v2["group"]]["treated"].append(s)
    else:
        if isinstance(v2["group"], list):
            for g in v2["group"]:
                GROUP2SAMPLE[g]["input"].append(s)
        else:
            GROUP2SAMPLE[v2["group"]]["input"].append(s)
    SAMPLE2BARCODE[s] = parse_barcode(v2.get("barcode", config["barcode"]))
    for i, v3 in enumerate(v2.get("data", []), 1):
        r = f"run{i}"
        SAMPLE2RUN[s][r] = {
            k4: os.path.expanduser(v4)
            if os.path.isabs(os.path.expanduser(v4))
            else os.path.relpath(os.path.expanduser(v4), WORKDIR)
            for k4, v4 in v3.items()
        }
        if not v2.get("forward_stranded", config["forward_stranded"]):
            SAMPLE2RUN[s][r]["R1"], SAMPLE2RUN[s][r]["R2"] = (
                SAMPLE2RUN[s][r]["R2"],
                SAMPLE2RUN[s][r]["R1"],
            )
    for k, v3 in v2.get("bam", {}).items():
        SAMPLE2BAM[s][k] = (
            os.path.expanduser(v3)
            if os.path.isabs(os.path.expanduser(v3))
            else os.path.relpath(os.path.expanduser(v3), WORKDIR)
        )


rule all:
    input:
        "report_reads/readsStats.html" if len(SAMPLE2RUN) > 0 else [],
        expand("call_sites/{reftype}.tsv.gz", reftype=REFTYPE),
        expand("filter_sites/{reftype}.tsv.gz", reftype=REFTYPE),


#### process reads ####


rule join_pairend_reads:
    input:
        lambda wildcards: SAMPLE2RUN[wildcards.sample][wildcards.rn].values(),
    output:
        fq=temp(os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}.fq.gz")),
        html="report_reads/joining/{sample}_{rn}.fastp.html",
        json="report_reads/joining/{sample}_{rn}.fastp.json",
    params:
        m=os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}_merge.fq.gz"),
        u1=os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}_u1.fq.gz"),
        u2=os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}_u2.fq.gz"),
        path_fastp=config["path"]["fastp"],
        path_joinFastq=config["path"]["joinFastq"],
    threads: 10
    run:
        if len(input) == 2:
            shell(
                """
        {params.path_fastp} --thread {threads} \
            --disable_adapter_trimming --merge --correction --overlap_len_require 10 --overlap_diff_percent_limit 20 \
            -i {input[0]} -I {input[1]} --merged_out {params.m} --out1 {params.u1} --out2 {params.u2} -h {output.html} -j {output.json}
        {params.path_joinFastq} {params.m} {params.u1} {params.u2} {output.fq}
        rm -f {params.m} {params.u1} {params.u2}
        """
            )
        else:
            shell(
                """
        ln -sfr {input[0]} {output.fq}
        touch {output.html} {output.json}
        """
            )


rule run_cutadapt:
    input:
        os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}.fq.gz"),
    output:
        fastq_trimmed=temp(
            os.path.join(TEMPDIR, "trimmed_reads/{sample}_{rn}_cut.fq.gz")
        ),
        fastq_short="discarded_reads/{sample}_{rn}_short.fq.gz"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_short.fq.gz"),
        report="report_reads/trimming/{sample}_{rn}_cutadapt.report",
    params:
        path_cutadapt=config["path"]["cutadapt"],
        p7=lambda wildcards: SAMPLE2BARCODE[wildcards.sample]["inline"]
        + config["adapter"]["p7"],
        drop_untrimmed_args=lambda wildcards: (
            "--untrimmed-output="
            + f"discarded_reads/{wildcards.sample}_{wildcards.rn}_untrimmed.fq.gz"
            if config["keep_discarded"]
            else "--discard-untrimmed"
        )
        if config["trimmed_only"] is True
        or (
            config["trimmed_only"] == "auto"
            and len(SAMPLE2BARCODE[wildcards.sample]["inline"]) > 0
        )
        else "",
        trim_p5_step=lambda wildcards, threads: " ".join(
            [
                config["path"]["cutadapt"],
                "-j",
                        str(threads),
                        "-g",
                        '"{};o=3;e=0.2;rightmost"'.format(config["adapter"]["p5"][-13:]),
                        "-",
                        "|",
                    ]
                )
                if config["trim_p5"]
                else "",
        trim_polyA_step=lambda wildcards, threads: " ".join(
            [
                config["path"]["cutadapt"],
                "-j",
                        str(threads),
                        "-a",
                        '"A{20};o=6;e=0.15"',
                        "-",
                        "|",
                    ]
                )
                if config["trim_polyA"]
                else "",
        extract_umi_args=lambda wildcards: "-u {} -u -{} ".format(
            SAMPLE2BARCODE[wildcards.sample]["umi5"],
            SAMPLE2BARCODE[wildcards.sample]["umi3"],
        )
        + ' --rename="{id}_{cut_prefix}{cut_suffix}"'
        if SAMPLE2BARCODE[wildcards.sample]["umi5"] > 0
        and SAMPLE2BARCODE[wildcards.sample]["umi3"] > 0
        else "-u {}".format(SAMPLE2BARCODE[wildcards.sample]["umi5"])
        + ' --rename="{id}_{cut_prefix}"'
        if SAMPLE2BARCODE[wildcards.sample]["umi5"] > 0
        else "-u -{}".format(SAMPLE2BARCODE[wildcards.sample]["umi3"])
        + ' --rename="{id}_{cut_suffix}"'
        if SAMPLE2BARCODE[wildcards.sample]["umi3"] > 0
        else "",
        mask_ends_args=lambda wildcards: "-u {}".format(
            SAMPLE2BARCODE[wildcards.sample]["mask5"]
        )
        if SAMPLE2BARCODE[wildcards.sample]["mask5"] > 0
        else "" + "-u -{}".format(SAMPLE2BARCODE[wildcards.sample]["mask3"])
        if SAMPLE2BARCODE[wildcards.sample]["mask3"] > 0
        else "",
    threads: 20
    shell:
        """
        {params.path_cutadapt} -j {threads} \
            --strip-suffix "/1" --strip-suffix "/2" \
            --strip-suffix ".1" --strip-suffix ".2" \
            -a "{params.p7};o=3;e=0.15" \
            {params.drop_untrimmed_args} \
            {input} | \
        {params.trim_p5_step} \
        {params.path_cutadapt} -j {threads} \
            {params.extract_umi_args} \
            - | \
        {params.trim_polyA_step} \
        {params.path_cutadapt} -j {threads} \
            {params.mask_ends_args} \
            -q 20 \
            --nextseq-trim=20  \
            --max-n=0 \
            -m 18 \
            --too-short-output {output.fastq_short} \
            - \
            -o {output.fastq_trimmed} \
            >{output.report}
        """


rule build_bowtie2_index:
    input:
        fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    output:
        idx=os.path.join(INTERNALDIR, "mapping_index/{reftype}.1.bt2")
        if config["keep_internal"]
        else temp(os.path.join(INTERNALDIR, "mapping_index/{reftype}.1.bt2")),
    params:
        path_bowtie2build=config["path"]["bowtie2Build"],
        ref_bowtie2=os.path.join(INTERNALDIR, "mapping_index/{reftype}"),
    threads: 2
    shell:
        """
        export LC_ALL=C
        {params.path_bowtie2build} --threads {threads} {input.fa} {params.ref_bowtie2}
        """


rule map_to_contamination_by_bowtie2:
    input:
        fq=os.path.join(TEMPDIR, "trimmed_reads/{sample}_{rn}_cut.fq.gz"),
        idx=lambda wildcards: REF["contamination"].get(
            "bt2", os.path.join(INTERNALDIR, "mapping_index/contamination")
        )
            + ".1.bt2",
    output:
        bam=temp(
            os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_contamination.bam")
        ),
        un=temp(os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_contamination.fq")),
        report="report_reads/mapping/{sample}_{rn}_contamination.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        path_samtools=config["path"]["samtools"],
        ref_bowtie2=lambda wildcards: REF["contamination"].get(
            "bt2", os.path.join(INTERNALDIR, "mapping_index/contamination")
        ),
        args_bowtie2="--sensitive -L 16 --rdg 0,2"
        if config["speedy_mapping"]
        else "--local --ma 2 --score-min G,20,8 -D 20 -R 3 -L 16 -N 1 --mp 4 --rdg 0,2"
        if config["greedy_mapping"]
        else "--end-to-end --ma 0 --score-min L,2,-0.5 -D 20 -R 3 -L 16 -N 1 --mp 4 --rdg 0,2",
    threads: 24
    shell:
        """
        export LC_ALL=C
        {params.path_bowtie2} -p {threads} \
            {params.args_bowtie2} \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {input.fq} 2>{output.report} | \
            {params.path_samtools} view -O BAM -o {output.bam}
        """


rule extract_contamination_unmap:
    input:
        os.path.join(TEMPDIR, "mapping_discarded/{sample}_{rn}_contamination.cram"),
    output:
        temp(os.path.join(TEMPDIR, "mapping_rerun/{sample}_{rn}_contamination.fq")),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=lambda wildcards: REF.get("contamination", {"fa": []})["fa"],
    threads: 4
    shell:
        """
        {params.path_samtools} fastq -@ {threads} --reference {params.ref_fa} {input} > {output}
        """


rule map_to_genes_by_bowtie2:
    input:
        fq=[
            os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_contamination.fq"),
            os.path.join(TEMPDIR, "mapping_rerun/{sample}_{rn}_contamination.fq"),
        ]
        if "contamination" in REF
        else os.path.join(TEMPDIR, "trimmed_reads/{sample}_{rn}_cut.fq.gz"),
        idx=lambda wildcards: REF["genes"].get(
            "bt2", os.path.join(INTERNALDIR, "mapping_index/genes")
        )
            + ".1.bt2",
    output:
        bam=temp(os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_genes.bam")),
        un=temp(os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_genes.fq")),
        report="report_reads/mapping/{sample}_{rn}_genes.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        path_samfilter=config["path"]["samfilter"],
        path_samtools=config["path"]["samtools"],
        ref_bowtie2=lambda wildcards: REF["genes"].get(
            "bt2", os.path.join(INTERNALDIR, "mapping_index/genes")
        ),
        args_bowtie2="--local --ma 2 --score-min G,10,7 -D 20 -R 3 -L 8 -N 1 -i S,1,0.5 --mp 6,3 --rdg 1,2 --rfg 6,3"
        if config["greedy_mapping"]
        else "--end-to-end --ma 0 --score-min L,4,-0.5 -D 20 -R 3 -L 8 -N 1 -i S,1,0.5 --mp 6,3 --rdg 1,2 --rfg 6,3",
        fq=lambda wildcards: os.path.join(
            TEMPDIR,
            f"mapping_unsort/{wildcards.sample}_{wildcards.rn}_contamination.fq",
        )
        + ","
        + os.path.join(
            TEMPDIR,
            f"mapping_rerun/{wildcards.sample}_{wildcards.rn}_contamination.fq",
        )
        if "contamination" in REF
        else os.path.join(
            TEMPDIR, f"trimmed_reads/{wildcards.sample}_{wildcards.rn}_cut.fq.gz"
        ),
    threads: 24
    shell:
        """
        export LC_ALL=C
        {params.path_bowtie2} -p {threads} \
            {params.args_bowtie2} --norc -a \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {params.fq} 2>{output.report} | \
            {params.path_samfilter} | \
            {params.path_samtools} view -O BAM -o {output.bam}
        """


rule extract_genes_unmap:
    input:
        os.path.join(TEMPDIR, "mapping_discarded/{sample}_{rn}_genes.cram"),
    output:
        temp(os.path.join(TEMPDIR, "mapping_rerun/{sample}_{rn}_genes.fq")),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=REF["genes"]["fa"],
    threads: 4
    shell:
        """
        {params.path_samtools} fastq -@ {threads} --reference {params.ref_fa} {input} > {output}
        """


rule map_to_genome_by_star:
    input:
        f1=os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_genes.fq"),
        f2=os.path.join(TEMPDIR, "mapping_rerun/{sample}_{rn}_genes.fq"),
    output:
        bam=temp(os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_genome.bam")),
        un="discarded_reads/{sample}_{rn}_unmapped.fq.gz"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_unmapped.fq.gz"),
        report="report_reads/mapping/{sample}_{rn}_genome.report",
        log_out=temp(os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_Log.out")),
        SJ_out=temp(os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_SJ.out.tab")),
        progress_out=temp(
            os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_Log.progress.out")
        ),
        std_out=temp(os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_Log.std.out")),
    params:
        output_pre=os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_"),
        un=os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_Unmapped.out.mate1"),
        report=os.path.join(TEMPDIR, "star_mapping/{sample}_{rn}_Log.final.out"),
        path_star=config["path"]["star"],
        ref_star=REF["genome"]["star"],
        match_prop=config["cutoff"]["min_match_prop"],
    threads: 24
    shell:
        """
        ulimit -n 20000
        rm -f {params.un}
        mkfifo {params.un}
        cat {params.un} | gzip > {output.un} &
        {params.path_star} \
          --runThreadN {threads} \
          --genomeDir {params.ref_star} \
          --readFilesIn {input.f1},{input.f2} \
          --alignEndsType Local \
          --scoreDelOpen -1 \
          --scoreDelBase -1 \
          --scoreInsOpen -2 \
          --scoreInsBase -2 \
          --outFilterMatchNmin 15 \
          --outFilterMatchNminOverLread {params.match_prop} \
          --outFilterMismatchNmax 10 \
          --outFilterMismatchNoverLmax 0.2 \
          --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
          --alignSJDBoverhangMin 1 \
          --alignSJoverhangMin 5 \
          --chimSegmentMin 20 \
          --chimOutType WithinBAM HardClip \
          --chimJunctionOverhangMin 15 \
          --chimScoreJunctionNonGTAG 0 \
          --outFilterMultimapNmax 10 \
          --outFilterMultimapScoreRange 0 \
          --outSAMmultNmax -1 \
          --outMultimapperOrder Random \
          --outReadsUnmapped Fastx \
          --outSAMtype BAM Unsorted \
          --outStd BAM_Unsorted \
          --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} LB:RNA PL:Illumina PU:SE \
          --outSAMattributes NH HI AS nM NM MD jM jI MC ch \
          --outFileNamePrefix {params.output_pre} > {output.bam}
        mv {params.report} {output.report}
        rm {params.un}
        """


rule gap_realign:
    input:
        os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_{reftype}.bam"),
    output:
        temp(
            os.path.join(
                TEMPDIR, "mapping_realigned_unsorted/{sample}_{rn}_{reftype}.cram"
            )
        ),
    params:
        path_realignGap=config["path"]["realignGap"],
        ref_fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    shell:
        """
        {params.path_realignGap} -r {params.ref_fa} -i {input} -o {output}
        """


rule sort_cal_filter_bam:
    input:
        os.path.join(TEMPDIR, "mapping_realigned_unsorted/{sample}_{rn}_{reftype}.cram"),
    output:
        cram=os.path.join(
            INTERNALDIR, "mapping_realigned/{sample}_{rn}_{reftype}.cram"
        )
        if config["keep_internal"]
        else temp(
            os.path.join(INTERNALDIR, "mapping_realigned/{sample}_{rn}_{reftype}.cram")
        ),
        un=temp(os.path.join(TEMPDIR, "mapping_discarded/{sample}_{rn}_{reftype}.cram")),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    threads: 8
    shell:
        """
        {params.path_samtools} sort -@ {threads} -m 4G {input} | \
            {params.path_samtools} calmd -@ {threads} - {params.ref_fa} 2>/dev/null | \
            {params.path_samtools} view -@ {threads} --reference {params.ref_fa} -e '[NM]<=5 && [NM]/(qlen-sclen)<=0.1' -O CRAM -U {output.un} -o {output.cram}
        """


rule combine_mapping_discarded:
    input:
        os.path.join(TEMPDIR, "mapping_discarded/{sample}_{rn}_genome.cram"),
    output:
        "discarded_reads/{sample}_{rn}_filteredmap.fq.gz"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_filteredmap.fq.gz"),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=REF["genome"]["fa"],
    threads: 4
    shell:
        """
        {params.path_samtools} fastq -@ {threads} --reference {params.ref_fa} -0 {output} {input}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            os.path.join(
                INTERNALDIR,
                f"mapping_realigned/{wildcards.sample}_{r}_{wildcards.reftype}.cram",
            )
            for r in SAMPLE2RUN[wildcards.sample]
        ],
    output:
        bam=temp(os.path.join(TEMPDIR, "combined_mapping/{sample}_{reftype}.bam")),
        bai=temp(os.path.join(TEMPDIR, "combined_mapping/{sample}_{reftype}.bam.bai")),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    threads: 8
    run:
        if len(input) > 1:
            shell(
                "{params.path_samtools} merge -@ {threads} --reference {params.ref_fa} --write-index -O BAM -o {output.bam}##idx##{output.bai} {input}"
            )
        else:
            shell(
                "{params.path_samtools} view -@ {threads} --reference {params.ref_fa} --write-index -O BAM -o {output.bam}##idx##{output.bai} {input}"
            )


rule drop_duplicates:
    input:
        bam=os.path.join(TEMPDIR, "combined_mapping/{sample}_{reftype}.bam"),
        bai=os.path.join(TEMPDIR, "combined_mapping/{sample}_{reftype}.bam.bai"),
    output:
        bam="align_bam/{sample}_{reftype}.bam",
        log="report_reads/deduping/{sample}_{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
        TEMPDIR=TEMPDIR,
    threads: 8
    run:
        if (
            SAMPLE2BARCODE[wildcards.sample]["umi5"]
            + SAMPLE2BARCODE[wildcards.sample]["umi3"]
            > 0
        ):
            shell(
                """
                java -server -Xmx46G -Xms24G -Xss100M -Djava.io.tmpdir={params.TEMPDIR} -jar {params.path_umicollapse} bam \
                    -t {threads} --data naive --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.log}
                """
            )
        else:
            shell(
                """
                cp {input.bam} {output.bam}
                touch {output.log}
                """
            )


rule index_dedup_bam:
    input:
        "align_bam/{sample}_{reftype}.bam",
    output:
        "align_bam/{sample}_{reftype}.bam.bai",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} index -@ {threads} {input}"


rule stat_dedup_bam:
    input:
        "align_bam/{sample}_{reftype}.bam",
    output:
        "report_reads/deduping/{sample}_{reftype}_dedup.report",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} flagstat -@ {threads} -O tsv {input} > {output}"


rule report_reads_stat:
    input:
        lambda wildcards: [
            f"report_reads/trimming/{s}_{r}_cutadapt.report"
            for s, v in SAMPLE2RUN.items()
            for r in v
        ],
        lambda wildcards: [
            f"report_reads/mapping/{s}_{r}_{t}.report"
            for s, v in SAMPLE2RUN.items()
            for r in v
            for t in REF.keys()
        ],
        lambda wildcards: [
            f"report_reads/deduping/{s}_{t}_dedup.report"
            for s in SAMPLE2RUN
            for t in REF.keys()
        ],
        # TODO: add discarded reads
        lambda wildcards: [
            f"discarded_reads/{s}_{r}_filteredmap.fq.gz"
            for s, v in SAMPLE2RUN.items()
            for r in v
        ],
    output:
        "report_reads/readsStats.html",
    params:
        path_multiqc=config["path"]["multiqc"],
    shell:
        "{params.path_multiqc} -f -m readsStats -t yc --no-data-dir -n {output} {input}"


##### call pU sites #####


rule merge_treated_bam_by_group:
    input:
        bam=lambda wildcards: [
            f"align_bam/{s}_{{reftype}}.bam"
            if s in SAMPLE2RUN
            else SAMPLE2BAM[s][wildcards.reftype]
            for s in GROUP2SAMPLE[wildcards.group]["treated"]
        ],
        bai=lambda wildcards: [
            f"align_bam/{s}_{{reftype}}.bam.bai"
            if s in SAMPLE2RUN
            else SAMPLE2BAM[s][wildcards.reftype] + ".bai"
            for s in GROUP2SAMPLE[wildcards.group]["treated"]
        ],
    output:
        bam=temp(os.path.join(TEMPDIR, "drop_duplicates_grouped/{group}_{reftype}.bam")),
        bai=temp(
            os.path.join(TEMPDIR, "drop_duplicates_grouped/{group}_{reftype}.bam.bai")
        ),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    shell:
        "{params.path_samtools} merge -@ {threads} --write-index -O BAM -o {output.bam}##idx##{output.bai} {input.bam}"


rule perbase_count_pre:
    input:
        bam=os.path.join(TEMPDIR, "drop_duplicates_grouped/{group}_{reftype}.bam"),
        bai=os.path.join(TEMPDIR, "drop_duplicates_grouped/{group}_{reftype}.bam.bai"),
    output:
        temp(os.path.join(TEMPDIR, "selected_region_by_group/{group}_{reftype}.bed")),
    params:
        path_delfilter=config["path"]["delfilter"],
        min_group_gap=config["cutoff"]["min_group_gap"],
        min_group_depth=config["cutoff"]["min_group_depth"],
        min_group_ratio=config["cutoff"]["min_group_ratio"],
    threads: 1
    shell:
        """
        {params.path_delfilter} -i {input.bam} -g {params.min_group_gap} -d {params.min_group_depth} -r {params.min_group_ratio} > {output}
        """


rule generate_faidx:
    input:
        fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    output:
        fai=os.path.join(INTERNALDIR, "fa_index/{reftype}.fa.fai")
        if config["keep_internal"]
        else temp(os.path.join(INTERNALDIR, "fa_index/{reftype}.fa.fai")),
    params:
        path_samtools=config["path"]["samtools"],
    shell:
        """
        {params.path_samtools} faidx {input.fa} --fai-idx {output.fai}
        """


rule prepare_bed_file:
    input:
        bed=expand(
            os.path.join(TEMPDIR, "selected_region_by_group/{group}_{{reftype}}.bed"),
            group=[g for g, s in GROUP2SAMPLE.items() if "treated" in s],
        ),
        fai=os.path.join(INTERNALDIR, "fa_index/{reftype}.fa.fai"),
    output:
        tmp=temp(os.path.join(TEMPDIR, "selected_region/picked_{reftype}_tmp.bed")),
        fwd=temp(os.path.join(TEMPDIR, "selected_region/picked_{reftype}_fwd.bed")),
        rev=temp(os.path.join(TEMPDIR, "selected_region/picked_{reftype}_rev.bed")),
    params:
        path_bedtools=config["path"]["bedtools"],
        min_group_num=config["cutoff"]["min_group_num"],
    threads: 4
    shell:
        """
        cat {input.bed} | {params.path_bedtools} slop -i - -g {input.fai} -b 3 | sort -S 4G --parallel={threads} -k1,1 -k2,2n >{output.tmp}
        {params.path_bedtools} merge -s -S + -c 1 -o count -i {output.tmp} | awk '$4 >= {params.min_group_num}' > {output.fwd}
        {params.path_bedtools} merge -s -S - -c 1 -o count -i {output.tmp} | awk '$4 >= {params.min_group_num}' > {output.rev}
        """


rule count_base_by_sample:
    input:
        bed=lambda wildcards: os.path.join(
            TEMPDIR,
            f"selected_region/picked_{wildcards.reftype}_{wildcards.orientation}.bed",
        )
        if wildcards.reftype in config["select_region"]
        else [],
        bam=lambda wildcards: "align_bam/{sample}_{reftype}.bam"
        if wildcards.sample in SAMPLE2RUN
        else SAMPLE2BAM[wildcards.sample][wildcards.reftype],
        bai=lambda wildcards: "align_bam/{sample}_{reftype}.bam.bai"
        if wildcards.sample in SAMPLE2RUN
        else SAMPLE2BAM[wildcards.sample][wildcards.reftype] + ".bai",
    output:
        temp(
            os.path.join(
                TEMPDIR, "pileup_bases_by_sample/{sample}_{reftype}_{orientation}.tsv"
            )
        ),
    params:
        path_samtools=config["path"]["samtools"],
        path_cpup=config["path"]["cpup"],
        ref=lambda wildcards: REF[wildcards.reftype]["fa"],
        region=lambda wildcards: "-l "
        + os.path.join(
            TEMPDIR,
            f"selected_region/picked_{wildcards.reftype}_{wildcards.orientation}.bed",
        )
        if wildcards.reftype in config["select_region"]
        else "",
        strand=lambda wildcards: "+" if wildcards.orientation == "fwd" else "-",
        flag=lambda wildcards: "--ff 3608"
        if wildcards.orientation == "fwd"
        else "--rf 16 --ff 3592",
    threads: 2
    shell:
        """
        {params.path_samtools} mpileup -aa -B -d 0 {params.flag} -Q 5 --reverse-del {params.region} -f {params.ref} {input.bam} | \
            {params.path_cpup} -H -S -i | \
            sed 's/\\t/\\t{params.strand}\\t/3' > {output}
        """


rule count_bases_combined:
    input:
        fwd=expand(
            os.path.join(
                TEMPDIR, "pileup_bases_by_sample/{sample}_{{reftype}}_fwd.tsv"
            ),
            sample=SAMPLE_IDS,
        ),
        rev=expand(
            os.path.join(
                TEMPDIR, "pileup_bases_by_sample/{sample}_{{reftype}}_rev.tsv"
            ),
            sample=SAMPLE_IDS,
        ),
    output:
        temp(os.path.join(TEMPDIR, "pileup_bases/{reftype}.tsv")),
    params:
        header="\t".join(["chr", "pos", "ref_base", "strand"] + list(SAMPLE_IDS)),
        idx=",".join(
            ["1", "2", "3", "4"] + [str(5 + i * 5) for i in range(len(SAMPLE_IDS))]
        ),
    threads: 4
    shell:
        """
        (
          echo {params.header:q}
          paste {input.fwd} | cut -f {params.idx}
          paste {input.rev} | cut -f {params.idx}
        ) >{output}
        """


rule adjust_sites:
    input:
        os.path.join(TEMPDIR, "pileup_bases/{reftype}.tsv"),
    output:
        "call_sites/{reftype}.tsv.gz",
    params:
        path_adjustGap=config["path"]["adjustGap"],
    shell:
        """
        {params.path_adjustGap} -i {input} -o {output}
        """


rule pre_filter_sites:
    input:
        "call_sites/{reftype}.tsv.gz",
    output:
        temp(os.path.join(TEMPDIR, "prefilter_sites/{reftype}.tsv.gz")),
    params:
        path_filterGap=config["path"]["filterGap"],
        min_group_gap=config["cutoff"]["min_group_gap"],
        min_group_depth=config["cutoff"]["min_group_depth"],
        min_group_ratio=config["cutoff"]["min_group_ratio"],
        min_group_num=config["cutoff"]["min_group_num"],
        columns=" ".join(
            [
                "-c "
                + ",".join(
                    [str(i) for i, s in enumerate(SAMPLE_IDS) if s in v["treated"]]
                )
                for v in GROUP2SAMPLE.values()
                if "treated" in v
            ]
        ),
    shell:
        """
        {params.path_filterGap} -i {input} -o {output} {params.columns} -g {params.min_group_gap} -d {params.min_group_depth} -r {params.min_group_ratio} -n {params.min_group_num}
        """


rule post_filter_sites:
    input:
        os.path.join(TEMPDIR, "prefilter_sites/{reftype}.tsv.gz"),
    output:
        "filter_sites/{reftype}.tsv.gz",
    params:
        group_filter=config.get("group_filter", {}),
        group_meta=dict(GROUP2SAMPLE),
        calibration_curve=CALI,
        ref_fasta=lambda wildcards: REF[wildcards.reftype]["fa"],
    script:
        config["path"]["pickSites"]
