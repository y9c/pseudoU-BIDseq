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


REF = config["reference"]
for k, v in REF.items():
    for k2, v2 in v.items():
        v2 = os.path.expanduser(v2)
        REF[k][k2] = v2 if os.path.isabs(v2) else os.path.relpath(v2, WORKDIR)

REFTYPE = ["genes", "genome"]
GROUP2SAMPLE = defaultdict(lambda: defaultdict(list))
SAMPLE_IDS = []
SAMPLE2RUN = defaultdict(dict)
SAMPLE2BAM = defaultdict(dict)
for s, v2 in config["samples"].items():
    SAMPLE_IDS.append(s)
    if v2.get("treated", True):
        GROUP2SAMPLE[v2["group"]]["treated"].append(s)
    else:
        GROUP2SAMPLE[v2["group"]]["input"].append(s)
    for i, v3 in enumerate(v2.get("data", []), 1):
        r = f"run{i}"
        SAMPLE2RUN[s][r] = {
            k4: os.path.relpath(os.path.expanduser(v4), WORKDIR)
            for k4, v4 in v3.items()
        }
    for k, v3 in v2.get("bam", {}).items():
        SAMPLE2BAM[s][k] = os.path.relpath(os.path.expanduser(v3), WORKDIR)


rule all:
    input:
        "report_reads/readsStats.html",
        expand("pileup_adjusted/{reftype}.tsv.gz", reftype=REFTYPE),


#### process reads ####


rule join_pairend_reads:
    input:
        lambda wildcards: SAMPLE2RUN[wildcards.sample][wildcards.rn].values(),
    output:
        temp(os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}.fq.gz")),
    params:
        path_fastp=config["path"]["fastp"],
        html="merged_reads/{sample}_{rn}.fastp.html",
        json="merged_reads/{sample}_{rn}.fastp.json",
    threads: 10
    run:
        if len(input) == 2:
            shell(
                """
        {params.path_fastp} --thread {threads} \
            --disable_adapter_trimming --merge --correction --overlap_len_require 10 --overlap_diff_percent_limit 20 \
            -i {input[0]} -I {input[1]} --merged_out {output} --out1 /dev/null --out2 /dev/null -h {params.html} -j {params.json}
        """
            )
        else:
            shell(
                """
        ln -sfr {input[0]} {output[0]}
        """
            )


rule run_cutadapt:
    input:
        os.path.join(TEMPDIR, "merged_reads/{sample}_{rn}.fq.gz"),
    output:
        fastq_trimmed=temp(
            os.path.join(TEMPDIR, "trimmed_reads/{sample}_{rn}_cut.fq.gz")
        ),
        fastq_untrimmed="discarded_reads/{sample}_{rn}_untrimmed.fq.gz"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_untrimmed.fq.gz"),
        fastq_short="discarded_reads/{sample}_{rn}_short.fq.gz"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_short.fq.gz"),
        report="report_reads/trimming/{sample}_{rn}_cutadapt.report",
    params:
        path_cutadapt=config["path"]["cutadapt"],
        p7=config["adapter"]["p7"],
        trim_p5_args='-n 3 -g "{};o=5;e=0.2" '.format(config["adapter"]["p5"])
        if config["trim_p5"]
        else "",
    threads: 20
    shell:
        """
        {params.path_cutadapt} -j {threads} \
            {params.trim_p5_args} \
            -a "{params.p7};o=3;e=0.15" \
            --untrimmed-output={output.fastq_untrimmed} \
            {input} | \
        {params.path_cutadapt} -j {threads} \
            -u -5 \
            --rename='{{id}}_{{cut_suffix}} {{comment}}' \
            -q 20 \
            --nextseq-trim=20  \
            --max-n=0 \
            -m 16 \
            --too-short-output {output.fastq_short} \
            - \
            -o {output.fastq_trimmed} \
            >{output.report}
        """


rule build_bowtie2_index:
    input:
        fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    output:
        idx=os.path.join(INTERNALDIR, "mapping_index/{reftype}.1.bt2"),
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
        args_bowtie2="--local --ma 2 --score-min G,20,8"
        if config["greedy_mapping"]
        else "--end-to-end --ma 0 --score-min L,2,-0.5",
    threads: 24
    shell:
        """
        export LC_ALL=C
        {params.path_bowtie2} -p {threads} \
            {params.args_bowtie2} -D 20 -R 3 -L 16 -N 1 --mp 4 --rdg 0,2 \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {input.fq} 2>{output.report} | \
            {params.path_samtools} view -O BAM -o {output.bam}
        """


rule map_to_genes_by_bowtie2:
    input:
        fq=os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_contamination.fq")
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
        args_bowtie2="--local --ma 2 --score-min G,10,7"
        if config["greedy_mapping"]
        else "--end-to-end --ma 0 --score-min L,4,-0.5",
    threads: 24
    shell:
        """
        export LC_ALL=C
        {params.path_bowtie2} -p {threads} \
            {params.args_bowtie2} --norc -D 20 -R 3 -L 8 -N 1 -i S,1,0.5 --mp 6,3 --rdg 0,2 -a \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {input.fq} 2>{output.report} | \
            {params.path_samfilter} | \
            {params.path_samtools} view -O BAM -o {output.bam}
        """


rule map_to_genome_by_star:
    input:
        os.path.join(TEMPDIR, "mapping_unsort/{sample}_{rn}_genes.fq"),
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
          --readFilesIn {input} \
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


rule fix_sort_filter_bam:
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
        un="discarded_reads/{sample}_{rn}_{reftype}_filteredmap.cram"
        if config["keep_discarded"]
        else temp("discarded_reads/{sample}_{rn}_{reftype}_filteredmap.cram"),
    params:
        path_samtools=config["path"]["samtools"],
        ref_fa=lambda wildcards: REF[wildcards.reftype]["fa"],
    threads: 8
    shell:
        """
        {params.path_samtools} calmd -@ {threads} {input} {params.ref_fa} | \
            {params.path_samtools} sort -@ {threads} -m 4G | \
            {params.path_samtools} view -@ {threads} --reference {params.ref_fa} -e '[NM]<=5' -O CRAM -U {output.un} -o {output.cram}
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
        bam="drop_duplicates/{sample}_{reftype}.bam",
        log="report_reads/deduping/{sample}_{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
        TEMPDIR=TEMPDIR,
    threads: 4
    shell:
        """
        java -server -Xmx46G -Xms24G -Xss100M -Djava.io.TEMPDIR={params.TEMPDIR} -jar {params.path_umicollapse} bam \
            --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.log}
        """


rule index_dedup_bam:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
    output:
        "drop_duplicates/{sample}_{reftype}.bam.bai",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} index -@ {threads} {input}"


rule stat_dedup_bam:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
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
            f"drop_duplicates/{s}_{{reftype}}.bam"
            if s in SAMPLE2RUN
            else SAMPLE2BAM[s][wildcards.reftype]
            for s in GROUP2SAMPLE[wildcards.group]["treated"]
        ],
        bai=lambda wildcards: [
            f"drop_duplicates/{s}_{{reftype}}.bam.bai"
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
        bed=os.path.join(TEMPDIR, "selected_region/picked_{reftype}_{orientation}.bed"),
        bam=lambda wildcards: "drop_duplicates/{sample}_{reftype}.bam"
        if wildcards.sample in SAMPLE2RUN
        else SAMPLE2BAM[wildcards.sample][wildcards.reftype],
        bai=lambda wildcards: "drop_duplicates/{sample}_{reftype}.bam.bai"
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
        strand=lambda wildcards: "+" if wildcards.orientation == "fwd" else "-",
        flag=lambda wildcards: "--ff 3608"
        if wildcards.orientation == "fwd"
        else "--rf 16 --ff 3592",
    threads: 2
    shell:
        """
        {params.path_samtools} mpileup -aa -B -d 0 {params.flag} -Q 5 --reverse-del -l {input.bed} -f {params.ref} {input.bam} | \
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
        temp(os.path.join(TEMPDIR, "pileup_bases/{reftype}.tsv.gz")),
    params:
        path_bgzip=config["path"]["bgzip"],
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
        ) | {params.path_bgzip} -@ {threads} -l 9 >{output}
        """


rule adjust_sites:
    input:
        os.path.join(TEMPDIR, "pileup_bases/{reftype}.tsv.gz"),
    output:
        "pileup_adjusted/{reftype}.tsv.gz",
    params:
        path_adjustGap=config["path"]["adjustGap"],
    shell:
        """
        {params.path_adjustGap} -i {input} -o {output}
        """
