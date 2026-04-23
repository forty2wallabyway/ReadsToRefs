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
FINAL_SCRIPT = "scripts/final_processing.py"

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
        "Figure1_EHDV_competition.png"

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
# Host removal
############################################

rule host_remove:
    input:
        reads=f"{WORK_DIR}/filtered/{{sample}}.fastq.gz",
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
# Adapter trimming
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
# Combine & index references (FIXED)
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
# Summary + heatmaps
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
# Final publication-quality figure (FIXED)
############################################

rule final_processing:
    input:
        summary=f"{REPORT_DIR}/alignment_summary_percentages.tsv"
    output:
        fig="Figure1_EHDV_competition.png"
    shell:
        """
        cd {REPORT_DIR}
        python ../../scripts/final_processing.py
        mv Figure1_EHDV_competition.png ../../{output.fig}
        """
