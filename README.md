## ReadsToRefs

### A workflow for competitive alignment of ONT reads against a set of predefined references.

#### This workflow was initially designed to analyze sequencing data collected during a collaborative project between CSU and the USDA; however, it can (in theory) be adapted to a variety of use cases.

#### Currently, it is designed to:
1. Perform quality filtering, remove host reads, subsample (optional), and trim adapters.
2. Competitively align reads against a set of provided references.
3. Output a handful of files, including: `alignment_heatmap_counts_top.png`, `alignment_heatmap_pct_top.png`, `alignment_summary_percentages.tsv`, and `parental_origin_annotated.tsv`

### The workflow assumes the following: 
1. A `data` directory containing raw fastq.gz reads.
2. A `references` directory with a `host.fasta` and any number (>2) of target `reference.fasta` files.
3. __Note:__ A `sispa_primers.fasta` file (located in `/references`) is needed for this particular workflow as there is a step that attempts to detect artifacts from library prep using SISPA + ONT.
   
#### Suggested usage:
1. Clone this repo to your desired workspace.
2. Create conda environment using provided `environment.yml` file.
3. Adjust desired presets in `config.yaml`
4. Activate workflow via `snakemake --cores [INT]`
