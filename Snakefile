# Snakefile
import os
import pandas as pd

SAMPLES_DF = pd.read_csv("data/samples.tsv", sep="\t")
SAMPLES = sorted(SAMPLES_DF["sample"].unique())

configfile: "config.yaml"

def get_fastq(sample, mate):
    """Get full fastq path for paired-end fastqs"""
    sub = SAMPLES_DF[SAMPLES_DF["sample"] == sample]
    # print(f"sample={sample}, mate={mate}, rows={len(sub)}")
    if mate == "R1":
        rows = sub[sub["filename"].str.contains("_R1")]
    else:
        rows = sub[sub["filename"].str.contains("_R2")]
    if len(rows) == 0:
        raise ValueError(f"no rows for sample={sample}, mate={mate}")
    row = rows.iloc[0]
    full = os.path.join(row["path"], row["filename"])
    # print(f"input({mate}) = {full}")
    return full

# def get_control(sample):
#     """Get control column from samples"""
#     sub = SAMPLES_DF[SAMPLES_DF["sample"] == sample]
#     ctrl_vals = sub["control"].dropna().unique()
#     if len(ctrl_vals) == 0:
#         return None
#     return ctrl_vals[0].strip() if ctrl_vals[0].strip() else None

rule all:
    input:
        expand("fastqFileQC/{sample}/{sample}_R1_fastqc.html", sample=SAMPLES),
        expand("fastqFileQC/{sample}/{sample}_R2_fastqc.html", sample=SAMPLES),
        expand("alignment/sam/{sample}_bowtie2.sam", sample=SAMPLES),
        expand("alignment/sam/{sample}_bowtie2_spikeIn.sam", sample=SAMPLES),
        expand("alignment/removeDuplicate/{sample}_bowtie2.sorted.rmDup.sam", sample=SAMPLES),
        expand("alignment/sam/fragmentLen/{sample}_fragmentLen.txt", sample=SAMPLES),
        expand("alignment/bam/{sample}_bowtie2.mapped.bam", sample=SAMPLES),
        expand("alignment/bed/{sample}_bowtie2.fragments.bed", sample=SAMPLES),
        expand("alignment/bed/{sample}_bowtie2.fragmentsCount.bin500.bed", sample=SAMPLES),
        expand("alignment/bedgraph/{sample}_bowtie2.fragments.normalized.bedgraph", sample=SAMPLES),
        expand("peakCalling/SEACR/{sample}_seacr_control.peaks.stringent.bed", sample=SAMPLES),
        expand("peakCalling/SEACR/{sample}_seacr_top0.01.peaks.stringent.bed", sample=SAMPLES)        

rule fastqc:
    input:
        R1=lambda w: get_fastq(w.sample, "R1"),
        R2=lambda w: get_fastq(w.sample, "R2")
    output:
        html_R1 = "fastqFileQC/{sample}/{sample}_R1_fastqc.html",
        zip_R1  = "fastqFileQC/{sample}/{sample}_R1_fastqc.zip",
        html_R2 = "fastqFileQC/{sample}/{sample}_R2_fastqc.html",
        zip_R2  = "fastqFileQC/{sample}/{sample}_R2_fastqc.zip"
    threads: config["threads"]["fastqc"]
    log:
        "logs/fastqc/{sample}.log"
    shell:
        r"""
        mkdir -p fastqFileQC/{wildcards.sample} logs/fastqc
        fastqc -t {threads} \
            -o fastqFileQC/{wildcards.sample} \
            "{input.R1}" "{input.R2}" &> {log}
        """

rule bowtie2_hg38:
    input:
        R1=lambda w: get_fastq(w.sample, "R1"),
        R2=lambda w: get_fastq(w.sample, "R2")
    output:
        sam="alignment/sam/{sample}_bowtie2.sam",
        summary="alignment/sam/bowtie2_summary/{sample}_bowtie2.txt"
    threads: config["threads"]["bowtie2"]
    params:
        ref=config["BOWTIE_INDEX"]
    log:
        "logs/bowtie2/{sample}_hg38.log"
    shell:
        r"""
        mkdir -p {config[alignment_dir]}/sam/bowtie2_summary logs/bowtie2
        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant \
            --phred33 -I 10 -X 700 -p {threads} \
            -x {params.ref} \
            -1 "{input.R1}" -2 "{input.R2}" \
            -S {output.sam} \
            &> {output.summary}
        """

rule bowtie2_spikein:
    input:
        R1=lambda w: get_fastq(w.sample, "R1"),
        R2=lambda w: get_fastq(w.sample, "R2")
    output:
        sam="alignment/sam/{sample}_bowtie2_spikeIn.sam",
        summary="alignment/sam/bowtie2_summary/{sample}_bowtie2_spikeIn.txt",
        seqdepth="alignment/sam/bowtie2_summary/{sample}_bowtie2_spikeIn.seqDepth"
    threads: config["threads"]["bowtie2"]
    params:
        ref=config["SPIKEIN_INDEX"]
    log:
        "logs/bowtie2/{sample}_spikeIn.log"
    shell:
        r"""
        mkdir -p {config[alignment_dir]}/sam/bowtie2_summary
        mkdir -p logs/bowtie2

        bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant \
            --phred33 -I 10 -X 700 -p {threads} \
            -x {params.ref} \
            -1 {input.R1} -2 {input.R2} \
            -S {output.sam} \
            &> {output.summary}

        # seqDepthDouble = number of aligned records (each pair counted twice),
        # seqDepth = paired-end fragments
        seqDepthDouble=$(samtools view -F 0x04 {output.sam} | wc -l)
        seqDepth=$((seqDepthDouble/2))
        echo $seqDepth > {output.seqdepth}
        """

rule picard_sortsam:
    input:
        sam="alignment/sam/{sample}_bowtie2.sam"
    output:
        sam_sorted="alignment/sam/{sample}_bowtie2.sorted.sam"
    params:
        picard_cmd=config["picard_cmd"]
    log:
        "logs/picard_sort/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/sam logs/picard_sort
        {params.picard_cmd} SortSam \
            I={input.sam} \
            O={output.sam_sorted} \
            SORT_ORDER=coordinate \
            TMP_DIR=./tmp \
            &> {log}
        """

rule picard_markdups:
    input:
        sam="alignment/sam/{sample}_bowtie2.sorted.sam"
    output:
        sam_dupmarked="alignment/removeDuplicate/{sample}_bowtie2.sorted.dupMarked.sam",
        metrics="alignment/removeDuplicate/picard_summary/{sample}_picard.dupMark.txt"
    params:
        picard_cmd=config["picard_cmd"]
    log:
        "logs/picard_markdups/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/removeDuplicate/picard_summary logs/picard_markdups
        {params.picard_cmd} MarkDuplicates \
            I={input.sam} \
            O={output.sam_dupmarked} \
            METRICS_FILE={output.metrics} \
            TMP_DIR=./tmp \
            &> {log}
        """

rule picard_remove_dupes:
    input:
        sam="alignment/sam/{sample}_bowtie2.sorted.sam"
    output:
        sam_rmdups="alignment/removeDuplicate/{sample}_bowtie2.sorted.rmDup.sam",
        metrics="alignment/removeDuplicate/picard_summary/{sample}_picard.rmDup.txt"
    params:
        picard_cmd=config["picard_cmd"]
    log:
        "logs/picard_rmdups/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/removeDuplicate/picard_summary logs/picard_rmdups
        {params.picard_cmd} MarkDuplicates \
            I={input.sam} \
            O={output.sam_rmdups} \
            METRICS_FILE={output.metrics} \
            REMOVE_DUPLICATES=true \
            TMP_DIR=./tmp \
            &> {log}
        """

rule fragment_len:
    input:
        sam="alignment/sam/{sample}_bowtie2.sam"
    output:
        fraglen="alignment/sam/fragmentLen/{sample}_fragmentLen.txt"
    log:
        "logs/fragment_len/{sample}.log"
    shell:
        r"""
        mkdir -p {config[alignment_dir]}/sam/fragmentLen logs/fragment_len
        samtools view -F 0x04 {input.sam} \
        | awk -F'\t' '
            function abs(x) {{ return ((x < 0.0) ? -x : x) }}
            {{ print abs($9) }}
        ' \
        | sort \
        | uniq -c \
        | awk -v OFS="\t" '{{print $2, $1/2}}' \
        > {output.fraglen} 2> {log}
        """

rule mapped_bam:
    input:
        sam="alignment/sam/{sample}_bowtie2.sam"
    output:
        bam="alignment/bam/{sample}_bowtie2.mapped.bam"
    threads: config["threads"]["bowtie2"]
    log:
        "logs/samtools_mapped/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/bam logs/samtools_mapped
        samtools view -bS -F 0x04 {input.sam} > {output.bam} 2> {log}
        """

rule bam_to_bed:
    input:
        bam="alignment/bam/{sample}_bowtie2.mapped.bam"
    output:
        bedpe = "alignment/bed/{sample}_bowtie2.bed",
        bed_clean = "alignment/bed/{sample}_bowtie2.clean.bed",
        fragments = "alignment/bed/{sample}_bowtie2.fragments.bed",
        bins_500bp = "alignment/bed/{sample}_bowtie2.fragmentsCount.bin500.bed"
    params:
        binLen = 500
    log:
        "logs/bam_to_bed/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/bed logs/bam_to_bed

        # Convert into bed file format
        bedtools bamtobed -i {input.bam} -bedpe > {output.bedpe} 2> {log}

        # Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
        awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {output.bedpe} > {output.bed_clean} 2>> {log}

        # Only extract the fragment related columns
        cut -f 1,2,6 {output.bed_clean} | sort -k1,1 -k2,2n -k3,3n > {output.fragments} 2>> {log}

        # 500 bp bins from fragment midpoints
        awk -v w={params.binLen} '{{print $1, int(($2 + $3)/(2*w))*w + w/2}}' {output.fragments} \
        | sort -k1,1V -k2,2n \
        | uniq -c \
        | awk -v OFS="\t" '{{print $2, $3, $1}}' \
        | sort -k1,1V -k2,2n \
        > {output.bins_500bp} 2>> {log}
        """

rule spikein_calibration:
    input:
        fragments="alignment/bed/{sample}_bowtie2.fragments.bed",
        seqdepth="alignment/sam/bowtie2_summary/{sample}_bowtie2_spikeIn.seqDepth"
    output:
        bedgraph="alignment/bedgraph/{sample}_bowtie2.fragments.normalized.bedgraph"
    params:
        chromSize=config["chromSize"]
    log:
        "logs/spikein_norm/{sample}.log"
    shell:
        r"""
        mkdir -p alignment/bedgraph logs/spikein_norm

        seqDepth=$(cat {input.seqdepth})
        if [ "$seqDepth" -gt 1 ]; then
            scale_factor=$(echo "10000 / $seqDepth" | bc -l)
            echo "Scaling factor for {wildcards.sample} is: $scale_factor" >> {log}

            bedtools genomecov -bg -scale "$scale_factor" -i {input.fragments} -g {params.chromSize} > {output.bedgraph} 2>> {log}
        else
            echo "seqDepth <= 1 for {wildcards.sample}, skipping normalization." >> {log}
            
            # create an empty file so Snakemake considers the target done
            touch {output.bedgraph}
        fi
        """

# Note next rule gave most trouble: for control (histControl) empty files created to bypass errors

rule seacr_peaks:
    input:
        bedgraph = "alignment/bedgraph/{sample}_bowtie2.fragments.normalized.bedgraph"
    output:
        peaks_ctrl = "peakCalling/SEACR/{sample}_seacr_control.peaks.stringent.bed",
        peaks_top = "peakCalling/SEACR/{sample}_seacr_top0.01.peaks.stringent.bed"
    params:
        seacr = config["seacr_script"],
        histControl = config.get("histControl", "").strip(),
        ctrl_suffix = "peakCalling/SEACR/{sample}_seacr_control.peaks",
        top_suffix = "peakCalling/SEACR/{sample}_seacr_top0.01.peaks"
    log:
        "logs/seacr/{sample}.log"
    shell:
        r"""
        mkdir -p peakCalling/SEACR logs/seacr

        tbg="{input.bedgraph}"
        seacr="{params.seacr}"
        histControl="{params.histControl}"

        # Always run top 0.01 mode
        bash "$seacr" "$tbg" 0.01 non stringent "{params.top_suffix}" >> {log} 2>&1

        # Run control mode if histControl != sample
        if [ -n "$histControl" ] && [ "$histControl" != "{wildcards.sample}" ]; then
            ctrl_bg="alignment/bedgraph/${{histControl}}_bowtie2.fragments.normalized.bedgraph"
            echo "Using control $histControl for {wildcards.sample} (bedgraph: $ctrl_bg)" >> {log}
            if [ -f "$ctrl_bg" ]; then
                bash "$seacr" "$tbg" "$ctrl_bg" non stringent "{params.ctrl_suffix}" >> {log} 2>&1
            else
                echo "Control bedgraph $ctrl_bg not found, creating empty control peaks file." >> {log}
                : > "{output.peaks_ctrl}"
            fi
        else
            echo "No valid control (or sample is its own control) for {wildcards.sample}, creating empty control peaks file." >> {log}
            : > "{output.peaks_ctrl}"
        fi
        """
