############################################
# Configuration
############################################

READS_DIR   = "data/reads"
REF_DIR     = "references"
ALIGN_DIR   = "alignments"
REPORTS_DIR = "reports"

############################################
# Auto-detect samples and references
############################################

SAMPLES, = glob_wildcards(f"{READS_DIR}/{{sample}}.fastq.gz")
REFS,    = glob_wildcards(f"{REF_DIR}/{{ref}}.fasta")

############################################
# Final target
############################################

rule all:
    input:
        f"{REPORTS_DIR}/alignment_summary_percentages.tsv",
        f"{REPORTS_DIR}/alignment_heatmap_top20.png"

############################################
# Index references for minimap2 (ONT)
############################################

rule index_reference:
    input:
        fasta = f"{REF_DIR}/{{ref}}.fasta"
    output:
        index = f"{REF_DIR}/{{ref}}.mmi"
    shell:
        """
        minimap2 -d {output.index} {input.fasta}
        """

############################################
# Align ONT reads with minimap2
############################################

rule align_reads_ont:
    input:
        reads = f"{READS_DIR}/{{sample}}.fastq.gz",
        ref   = f"{REF_DIR}/{{ref}}.mmi"
    output:
        bam = f"{ALIGN_DIR}/{{sample}}.{{ref}}.sorted.bam"
    threads: 8
    shell:
        """
        minimap2 \
            -x map-ont \
            -a \
            -t {threads} \
            {input.ref} \
            {input.reads} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam}
        """

############################################
# Count mapped ONT reads (primary only)
############################################

rule count_aligned_reads:
    input:
        bam = f"{ALIGN_DIR}/{{sample}}.{{ref}}.sorted.bam"
    output:
        txt = f"{REPORTS_DIR}/counts/{{sample}}.{{ref}}.txt"
    shell:
        """
        # Count primary, mapped reads
        samtools view \
            -c \
            -F 0x104 \
            {input.bam} \
            > {output.txt}
        """

############################################
# Per-sample percentage summary
############################################

rule summarize_counts_percentages:
    input:
        expand(
            f"{REPORTS_DIR}/counts/{{sample}}.{{ref}}.txt",
            sample=SAMPLES,
            ref=REFS
        )
    output:
        summary = f"{REPORTS_DIR}/alignment_summary_percentages.tsv"
    run:
        from collections import defaultdict

        counts = defaultdict(dict)
        totals = defaultdict(int)

        for sample in SAMPLES:
            for ref in REFS:
                fname = f"{REPORTS_DIR}/counts/{sample}.{ref}.txt"
                with open(fname) as f:
                    c = int(f.read().strip())
                counts[sample][ref] = c
                totals[sample] += c

        with open(output.summary, "w") as out:
            out.write("sample\treference\taligned_reads\tpercent_of_sample\n")
            for sample in SAMPLES:
                for ref in REFS:
                    aligned = counts[sample][ref]
                    pct = (aligned / totals[sample] * 100) if totals[sample] > 0 else 0
                    out.write(
                        f"{sample}\t{ref}\t{aligned}\t{pct:.2f}\n"
                    )

############################################
# Heatmap of top 20 references
############################################

rule alignment_heatmap:
    input:
        summary = f"{REPORTS_DIR}/alignment_summary_percentages.tsv"
    output:
        plot = f"{REPORTS_DIR}/alignment_heatmap_top20.png"
    shell:
        r"""
        python - <<'EOF'
        import pandas as pd
        import matplotlib.pyplot as plt
        import seaborn as sns

        df = pd.read_csv("{input.summary}", sep="\t")

        matrix = df.pivot(
            index="sample",
            columns="reference",
            values="percent_of_sample"
        ).fillna(0)

        top_refs = (
            matrix.mean(axis=0)
            .sort_values(ascending=False)
            .head(20)
            .index
        )

        matrix = matrix[top_refs]

        plt.figure(figsize=(1 + 0.5 * len(top_refs), 6))
        sns.heatmap(
            matrix,
            cmap="viridis",
            vmin=0,
            vmax=100,
            cbar_kws={"label": "Percent of aligned reads"},
            linewidths=0.5,
            linecolor="white"
        )

        plt.title("ONT read alignment heatmap (Top 20 references)")
        plt.xlabel("Reference")
        plt.ylabel("Sample")

        plt.tight_layout()
        plt.savefig("{output.plot}", dpi=300)
        plt.close()
        EOF
        """