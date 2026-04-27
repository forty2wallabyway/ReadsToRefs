from glob import glob

############################################
# Configuration
############################################

configfile: "config.yaml"

THREADS = config["threads"]

READS_DIR   = "data"
REFS_DIR    = "references"
WORK_DIR    = "work"
RESULTS_DIR = "results"
ALIGN_DIR   = f"{RESULTS_DIR}/alignments"
COUNT_DIR   = f"{RESULTS_DIR}/counts"
REPORT_DIR  = f"{RESULTS_DIR}/reports"

PLOT_SCRIPT  = "scripts/summarize_and_plot.py"

############################################
# Samples
############################################

SAMPLES, = glob_wildcards(f"{READS_DIR}/{{sample}}.fastq.gz")

if not SAMPLES:
    raise ValueError("No FASTQ files found in data/")

############################################
# References (EXCLUDE combined reference)
############################################

REF_FASTAS = glob(f"{REFS_DIR}/*.fasta")
REF_FASTAS = [f for f in REF_FASTAS if not f.endswith("all_refs.fasta")]

HOST_FASTA = f"{REFS_DIR}/host.fasta"
HOST_INDEX = f"{REFS_DIR}/host.mmi"

############################################
# Final targets
############################################

rule all:
    input:
        f"{REPORT_DIR}/alignment_summary_percentages.tsv",
        f"{REPORT_DIR}/alignment_heatmap_pct_top.png",
        f"{REPORT_DIR}/alignment_heatmap_counts_top.png",
        f"{RESULTS_DIR}/parental_origin.tsv",
        expand(f"{REPORT_DIR}/chimeras/{{sample}}.split_reads.tsv", sample=SAMPLES),
        expand(f"{REPORT_DIR}/chimeras/{{sample}}.internal_primers.tsv", sample=SAMPLES),
        f"{RESULTS_DIR}/parental_origin_annotated.tsv",
        expand(f"{RESULTS_DIR}/coverage/{{sample}}.coverage.tsv", sample=SAMPLES)

############################################
# Host reference indexing
############################################

rule index_host:
    input:
        HOST_FASTA
    output:
        HOST_INDEX
    shell:
        "minimap2 -d {output} {input}"

############################################
# Quality filtering
############################################

rule filter_reads:
    input:
        f"{READS_DIR}/{{sample}}.fastq.gz"
    output:
        temp(f"{WORK_DIR}/filtered/{{sample}}.fastq.gz")
    shell:
        """
        nanoq \
          --min-len {config[quality_filter][min_len]} \
          --min-qual {config[quality_filter][min_qual]} \
          -i {input} -o {output}
        """

############################################
## SISPA primer trimming
############################################

SISPA_PRIMERS = "references/sispa_primers.fasta"

rule trim_sispa_primers:
    input:
        f"{WORK_DIR}/filtered/{{sample}}.fastq.gz"
    output:
        temp(f"{WORK_DIR}/noprimer/{{sample}}.fastq.gz"),
        f"{REPORT_DIR}/sispa/{{sample}}.cutadapt.txt"
    threads: THREADS
    shell:
        """
        cutadapt \
          -g file:{SISPA_PRIMERS} \
          -a file:{SISPA_PRIMERS} \
          --discard-untrimmed \
          --report=full \
          -j {threads} \
          -o {output[0]} \
          {input} > {output[1]}
        """

############################################
# Host removal
############################################

rule host_remove:
    input:
        reads=f"{WORK_DIR}/noprimer/{{sample}}.fastq.gz",
        index=HOST_INDEX
    output:
        temp(f"{WORK_DIR}/nohost/{{sample}}.fastq.gz")
    threads: THREADS
    shell:
        """
        minimap2 -ax map-ont -t {threads} {input.index} {input.reads} |
        samtools view -b -f 4 |
        samtools fastq |
        gzip -c > {output}
        """

############################################
# Optional subsampling
############################################

rule subsample:
    input:
        f"{WORK_DIR}/nohost/{{sample}}.fastq.gz"
    output:
        temp(f"{WORK_DIR}/subsampled/{{sample}}.fastq.gz")
    run:
        if config["subsample"]["enabled"]:
            shell(
                "seqtk sample -s {seed} {input} {n} | gzip -c > {output}".format(
                    seed=config["subsample"]["seed"],
                    n=config["subsample"]["n_reads"],
                    input=input,
                    output=output,
                )
            )
        else:
            shell("cp {input} {output}")

############################################
# Barcode and ONT adapter trimming
############################################

rule trim_reads:
    input:
        f"{WORK_DIR}/subsampled/{{sample}}.fastq.gz"
    output:
        f"{RESULTS_DIR}/reads/{{sample}}.final.fastq.gz"
    threads: THREADS
    shell:
        "porechop -i {input} -o {output} --threads {threads}"

############################################
# Combine & index references
############################################

rule combine_refs:
    input:
        REF_FASTAS
    output:
        f"{REFS_DIR}/all_refs.fasta"
    shell:
        "cat {input} > {output}"

rule index_refs:
    input:
        f"{REFS_DIR}/all_refs.fasta"
    output:
        f"{REFS_DIR}/all_refs.mmi"
    shell:
        "minimap2 -d {output} {input}"

############################################
# Competitive alignment
############################################

rule align_competitive:
    input:
        reads=f"{RESULTS_DIR}/reads/{{sample}}.final.fastq.gz",
        ref=f"{REFS_DIR}/all_refs.mmi"
    output:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    threads: THREADS
    shell:
        """
        minimap2 -x {config[competitive][preset]} -a -t {threads} {input.ref} {input.reads} |
        samtools sort -@ {threads} -o {output}
        samtools index {output}
        """

############################################
# Calculate coverage
############################################

rule segment_coverage:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    output:
        f"{RESULTS_DIR}/coverage/{{sample}}.coverage.tsv"
    shell:
        """
        samtools coverage {input.bam} > {output}
        """

############################################
# Best-hit competitive counting
############################################

rule count_best_hit:
    input:
        f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    output:
        f"{COUNT_DIR}/{{sample}}.per_reference.tsv"
    shell:
        """
        samtools view -h -F 0x100 {input} |
        awk '
        BEGIN {{OFS="\\t"}}
        /^@/ {{next}}
        {{
          read=$1; ref=$3; mapq=$5; as=0
          for (i=12;i<=NF;i++)
            if ($i ~ /^AS:i:/) {{split($i,a,":"); as=a[3]}}
          if (!(read in best) || as>best_as[read] || (as==best_as[read] && mapq>best_mapq[read])) {{
            best[read]=ref; best_as[read]=as; best_mapq[read]=mapq
          }}
        }}
        END {{
          for (r in best)
            if (best[r]!="*") counts[best[r]]++
          for (ref in counts) print ref, counts[ref]
        }}' > {output}
        """

############################################
## Chimera detection
############################################

rule detect_split_alignments:
    input:
        bam=f"{ALIGN_DIR}/{{sample}}.sorted.bam"
    output:
        f"{REPORT_DIR}/chimeras/{{sample}}.split_reads.tsv"
    shell:
        """
        samtools view -h {input.bam} | \
        awk '
        /^@/ {{next}}
        {{
          if ($6 ~ /N/ || $6 ~ /[0-9]+S.*[0-9]+S/) {{
            print $1, $3, $6
          }}
        }}' OFS="\t" > {output}
        """
        
rule detect_internal_primers:
    input:
        f"{WORK_DIR}/noprimer/{{sample}}.fastq.gz"
    output:
        f"{REPORT_DIR}/chimeras/{{sample}}.internal_primers.tsv"
    shell:
        r"""
        cutadapt \
          -b file:{SISPA_PRIMERS} \
          --action=none \
          --error-rate 0.2 \
          --overlap 12 \
          --info-file {output}.info \
          -o /dev/null \
          {input}

        cut -f1 {output}.info | tail -n +2 | sort -u > {output}
        rm {output}.info
        """

############################################
# Summary and heatmaps
############################################

rule summarize_and_plot:
    input:
        expand(f"{COUNT_DIR}/{{sample}}.per_reference.tsv", sample=SAMPLES)
    output:
        summary=f"{REPORT_DIR}/alignment_summary_percentages.tsv",
        heatmap_pct=f"{REPORT_DIR}/alignment_heatmap_pct_top.png",
        heatmap_counts=f"{REPORT_DIR}/alignment_heatmap_counts_top.png"
    shell:
        """
        python {PLOT_SCRIPT} \
          --mode competitive \
          --top-n {config[competitive][top_n]} \
          --reads-dir {RESULTS_DIR}/reads \
          --reads-suffix .final.fastq.gz \
          --include-unmapped \
          --out-summary {output.summary} \
          --out-heatmap-pct {output.heatmap_pct} \
          --out-heatmap-counts {output.heatmap_counts} \
          {input}
        """

############################################
## Parental origin table
############################################

rule parental_origin:
    input:
        expand(f"{COUNT_DIR}/{{sample}}.per_reference.tsv", sample=SAMPLES)
    output:
        f"{RESULTS_DIR}/parental_origin.tsv"
    shell:
        """
        python scripts/parental_origin.py \
          --inputs {input} \
          --min-diff {config[parental][min_diff]} \
          --out {output}
        """

############################################
# Annotate chimeras
############################################

rule annotate_chimeras:
    input:
        parental=f"{RESULTS_DIR}/parental_origin.tsv",
        split=expand(f"{REPORT_DIR}/chimeras/{{sample}}.split_reads.tsv", sample=SAMPLES),
        primer=expand(f"{REPORT_DIR}/chimeras/{{sample}}.internal_primers.tsv", sample=SAMPLES)
    output:
        f"{RESULTS_DIR}/parental_origin_annotated.tsv"
    shell:
        """
        python scripts/annotate_chimeras.py \
          --parental {input.parental} \
          --bam-dir {ALIGN_DIR} \
          --split {input.split} \
          --primer {input.primer} \
          --out {output}
        """  
