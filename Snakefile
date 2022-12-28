import sys
from snakemake.utils import min_version
from collections import defaultdict

if sys.version_info < (3, 6):
    sys.exit("Python 3.6 or later is required.\n")
min_version("7.0")


configfile: "config.yaml"


BATCH = config["batch"]
SPECIES = config["species"]
WORKSPACE = config.get("workspace", f"workspace_{BATCH}")


workdir: WORKSPACE


REF = config["reference"]
REFTYPE = ["genes", "genome"]
REFTYPE_ALL = ["contamination", "genes", "genome"]
GROUP2SAMPLE = defaultdict(lambda: defaultdict(list))
SAMPLE_IDS = []
SAMPLE2RUN = defaultdict(list)
RUN2DATA = {}
SAMPLE2BAM = defaultdict(dict)
for n2, v2 in config[f"samples"].items():
    if v2.get("treated", True):
        GROUP2SAMPLE[v2["group"]]["treated"].append(n2)
    else:
        GROUP2SAMPLE[v2["group"]]["input"].append(n2)
    s = n2
    SAMPLE_IDS.append(s)
    for i, v3 in enumerate(v2.get("data", []), 1):
        r = f"{s}_run{i}"
        SAMPLE2RUN[s].append(r)
        RUN2DATA[r] = {k4: os.path.expanduser(v4) for k4, v4 in v3.items()}
    for k, v3 in v2.get("bam", {}).items():
        SAMPLE2BAM[s][k] = os.path.expanduser(v3)
META_DATA = {
    "batch": BATCH,
    "species": SPECIES,
    "ref": REF,
    "REFTYPE": REFTYPE_ALL,
    "sample_ids": list(SAMPLE2RUN.keys()),
    "run_ids": list(RUN2DATA.keys()),
    "groups": GROUP2SAMPLE,
}
META_FILE = "meta.json"
if not os.path.exists(META_FILE) or json.load(open(META_FILE, "r")) != META_DATA:
    with open(META_FILE, "w") as f:
        json.dump(META_DATA, f, indent=4)


rule all:
    input:
        expand("pileup_adjusted/{reftype}.tsv.gz", reftype=REFTYPE),


#### process reads ####


rule join_pairend_reads:
    input:
        lambda wildcards: RUN2DATA[wildcards.rn].values(),
    output:
        temp("merged_reads/{rn}.fq.gz"),
    params:
        path_fastp=config["path"]["fastp"],
        html="merged_reads/{rn}.fastp.html",
        json="merged_reads/{rn}.fastp.json",
    threads: 10
    run:
        if len(input) == 2:
            shell(
                """
        {params.path_fastp} --thread {threads} --disable_adapter_trimming --merge --correction --overlap_len_require 10 --overlap_diff_percent_limit 20 -i {input[0]} -I {input[1]} --merged_out {output} --out1 /dev/null --out2 /dev/null -h {params.html} -j {params.json}
        """
            )
        else:
            shell(
                """
        ln -sf {input[0]} {output}
        """
            )


rule cutadapt:
    input:
        "merged_reads/{rn}.fq.gz",
    output:
        fastq_trimmed="cut_adapter/{rn}_cut.fq.gz",
        fastq_untrimmed="cut_adapter/{rn}_untrimmed.fq.gz",
        fastq_short="cut_adapter/{rn}_short.fq.gz",
        report="cut_adapter/{rn}_cutadapt.report",
    params:
        path_cutadapt=config["path"]["cutadapt"],
        adapter3=config["adapter"]["p7"],
    threads: 20
    shell:
        """
        {params.path_cutadapt} -j {threads} \
            -a "{params.adapter3};o=3;e=0.15" \
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


rule map_to_contamination_by_bowtie2:
    input:
        "cut_adapter/{rn}_cut.fq.gz",
    output:
        sam=temp("mapping_unsort/{rn}_contamination.sam"),
        un=temp("mapping_unsort/{rn}_contamination.fq"),
        report="mapping_report/{rn}_contamination.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: REF["contamination"]["bt2"],
    threads: 24
    shell:
        """
        {params.path_bowtie2} -p {threads} \
            --end-to-end -D 20 -R 3 --score-min L,5,-0.5 -L 16 -N 1 --mp 4 --rdg 0,2 \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {input} > {output.sam} 2>{output.report}
        """


rule map_to_genes_by_bowtie2:
    input:
        "mapping_unsort/{rn}_contamination.fq",
    output:
        sam=temp("mapping_unsort/{rn}_genes.sam"),
        un=temp("mapping_unsort/{rn}_genes.fq"),
        report="mapping_report/{rn}_genes.report",
    params:
        path_bowtie2=config["path"]["bowtie2"],
        ref_bowtie2=lambda wildcards: REF["genes"]["bt2"],
        path_samfilter=config["path"]["samfilter"],
    threads: 24
    shell:
        """
        {params.path_bowtie2} -p {threads} \
            --end-to-end --norc -D 20 -R 3 --score-min L,4,-0.5 -L 10 -i S,1,0.5 -N 1 --mp 6,3 --rdg 0,2 -a \
            --no-unal --un {output.un} -x {params.ref_bowtie2} -U {input} 2>{output.report} | {params.path_samfilter} > {output.sam}
        """


rule map_to_genome_by_star:
    input:
        "mapping_unsort/{rn}_genes.fq",
    output:
        sam=temp("mapping_unsort/{rn}_genome.sam"),
        un="final_unmapped/{rn}.fq.gz",
        report="mapping_report/{rn}_genome.report",
        log="star_mapping/{rn}_Log.out",
    params:
        output_pre="star_mapping/{rn}_",
        sam="star_mapping/{rn}_Aligned.out.sam",
        un="star_mapping/{rn}_Unmapped.out.mate1",
        report="star_mapping/{rn}_Log.final.out",
        path_star=config["path"]["star"],
        ref_star=REF["genome"]["star"],
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
          --outFilterMatchNmin 15 \
          --outFilterMatchNminOverLread 0.66 \
          --outFilterMismatchNmax 10 \
          --outFilterMismatchNoverLmax 0.2 \
          --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
          --scoreDelOpen -1 \
          --scoreDelBase -1 \
          --scoreInsOpen -2 \
          --scoreInsBase -2 \
          --alignSJoverhangMin 5 \
          --outFilterMultimapNmax 10 \
          --outSAMmultNmax -1 \
          --outReadsUnmapped Fastx \
          --outSAMattrRGline ID:{wildcards.rn} SM:{wildcards.rn} LB:RNA PL:Illumina PU:SE \
          --outSAMattributes NH HI AS nM NM MD jM jI MC \
          --limitBAMsortRAM 8000000000 \
          --outSAMtype SAM \
          --outFileNamePrefix {params.output_pre}
        mv {params.report} {output.report}
        mv {params.sam} {output.sam}
        rm {params.un}
        """


rule mapping_sam_to_cram:
    input:
        "mapping_unsort/{rn}_{reftype}.sam",
    output:
        "mapping_sorted/{rn}_{reftype}.cram",
    params:
        ref_fa=REF["genome"]["fa"],
        path_samtools=config["path"]["samtools"],
    threads: 12
    shell:
        "{params.path_samtools} sort -@ {threads} --write-index --reference {params.ref_fa} -O CRAM -o {output} {input}"


rule gap_realign:
    input:
        "mapping_sorted/{rn}_{reftype}.cram",
    output:
        temp("mapping_realigned_tmp/{rn}_{reftype}.bam"),
    params:
        path_realignGap=config["path"]["realignGap"],
        ref_fa=(
            lambda wildcards: REF[wildcards.reftype]["fa"]
            if wildcards.reftype in REF
            else config[f"ref_{wildcards.reftype}"]["fa"]
        ),
    shell:
        """
        {params.path_realignGap} -r {params.ref_fa} -i {input} -o {output}
        """


rule sort_and_filter_bam:
    input:
        "mapping_realigned_tmp/{rn}_{reftype}.bam",
    output:
        cram="mapping_realigned/{rn}_{reftype}.cram",
        crai="mapping_realigned/{rn}_{reftype}.cram.crai",
    wildcard_constraints:
        reftype="contamination|genes|genome",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    shell:
        """
        {params.path_samtools} sort -@ {threads} --write-index --input-fmt-option 'filter=[NM]<=10' -m 4G -O CRAM -o {output.cram}##idx##{output.crai} {input}
        """


rule combine_runs:
    input:
        lambda wildcards: [
            os.path.join("mapping_realigned", f"{r}_{wildcards.reftype}.cram")
            for r in SAMPLE2RUN[wildcards.sample]
        ],
    output:
        bam=temp("combined_mapping/{sample}_{reftype}.bam"),
        bai=temp("combined_mapping/{sample}_{reftype}.bam.bai"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} merge -@ {threads} --write-index -O BAM -o {output.bam}##idx##{output.bai} {input}"


rule drop_duplicates:
    input:
        bam="combined_mapping/{sample}_{reftype}.bam",
        bai="combined_mapping/{sample}_{reftype}.bam.bai",
    output:
        bam="drop_duplicates/{sample}_{reftype}.bam",
        log="drop_duplicates/{sample}_{reftype}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
        tmpdir=config["tmp_dir"],
    threads: 4
    shell:
        """
        java -server -Xmx46G -Xms24G -Xss100M -Djava.io.tmpdir={params.tmpdir} -jar {params.path_umicollapse} bam \
            --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.log}
        """


rule index_bam_dedup:
    input:
        "drop_duplicates/{sample}_{reftype}.bam",
    output:
        "drop_duplicates/{sample}_{reftype}.bam.bai",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    shell:
        "{params.path_samtools} index -@ {threads} {input}"


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
        bam=temp("drop_duplicates_grouped/{group}_{reftype}.bam"),
        bai=temp("drop_duplicates_grouped/{group}_{reftype}.bam.bai"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    shell:
        "{params.path_samtools} merge -@ {threads} --write-index -O BAM -o {output.bam}##idx##{output.bai} {input.bam}"


rule perbase_count_pre:
    input:
        bam="drop_duplicates_grouped/{group}_{reftype}.bam",
        bai="drop_duplicates_grouped/{group}_{reftype}.bam.bai",
    output:
        temp("selected_region_by_group/{group}_{reftype}.bed"),
    params:
        path_delfilter=config["path"]["delfilter"],
    threads: 1
    shell:
        """
        {params.path_delfilter} -i {input.bam} -g 5 -d 10 -r 0.02 > {output}
        """


rule prepare_bed_file:
    input:
        expand(
            "selected_region_by_group/{group}_{{reftype}}.bed",
            group=[g for g, s in GROUP2SAMPLE.items() if "treated" in s],
        ),
    output:
        tmp=temp("selected_region/picked_{reftype}_tmp.bed"),
        fwd="selected_region/picked_{reftype}_fwd.bed",
        rev="selected_region/picked_{reftype}_rev.bed",
    params:
        fai=lambda wildcards: REF[wildcards.reftype]["fai"],
        path_bedtools=config["path"]["bedtools"],
        min_group=config["min_group"],
    threads: 4
    shell:
        """
        cat {input} | {params.path_bedtools} slop -i - -g {params.fai} -b 3 | sort -S 4G --parallel={threads} -k1,1 -k2,2n >{output.tmp}
        {params.path_bedtools} merge -s -S + -c 1 -o count -i {output.tmp} | awk '$4 >= {params.min_group}' > {output.fwd}
        {params.path_bedtools} merge -s -S - -c 1 -o count -i {output.tmp} | awk '$4 >= {params.min_group}' > {output.rev}
        """


rule count_base_by_sample:
    input:
        bed="selected_region/picked_{reftype}_{orientation}.bed",
        bam=lambda wildcards: "drop_duplicates/{sample}_{reftype}.bam"
        if wildcards.sample in SAMPLE2RUN
        else SAMPLE2BAM[wildcards.sample][wildcards.reftype],
        bai=lambda wildcards: "drop_duplicates/{sample}_{reftype}.bam.bai"
        if wildcards.sample in SAMPLE2RUN
        else SAMPLE2BAM[wildcards.sample][wildcards.reftype] + ".bai",
    output:
        "pileup_bases_by_sample/{sample}_{reftype}_{orientation}.tsv",
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
        {params.path_samtools} mpileup -aa -B -d 0 {params.flag} -Q 5 --reverse-del -l {input.bed} -f {params.ref} {input.bam} | {params.path_cpup} -H -S -i | sed 's/\\t/\\t{params.strand}\\t/3' > {output}
        """


rule count_bases_combined:
    input:
        fwd=expand(
            "pileup_bases_by_sample/{sample}_{{reftype}}_fwd.tsv", sample=SAMPLE_IDS
        ),
        rev=expand(
            "pileup_bases_by_sample/{sample}_{{reftype}}_rev.tsv", sample=SAMPLE_IDS
        ),
    output:
        "pileup_bases/{reftype}.tsv.gz",
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
        "pileup_bases/{reftype}.tsv.gz",
    output:
        "pileup_adjusted/{reftype}.tsv.gz",
    params:
        path_adjustGap=config["path"]["adjustGap"],
    shell:
        """
        {params.path_adjustGap} -i {input} -o {output}
        """
