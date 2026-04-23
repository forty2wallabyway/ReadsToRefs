## ReadsToRefs

### A workflow for competitive alignment of ONT (Nanopore) reads against a set of predefined references.

#### Currently designed to:
1. Perform quality filtering, remove host reads, subsample (optional), and trim adapters.
2. Competitively align reads against provided references (two minimum).
3. Output an `alignment_summary_percentages.tsv`, two heatmaps, and a final figure (in development).

### The workflow assumes the following directories are locally present: 
1. A `data` directory containing raw fastq.gz reads.
2. A `references` directory with a `host.fasta` and any number of target `reference.fasta` files.

#### Suggested usage:
1. Clone this repo to your desired workspace.
2. Create conda environment using provided environment file.
3. Adjust desired presets in `config.yaml`
4. Run workflow by activating Snakemake with the desired number of `--cores`
