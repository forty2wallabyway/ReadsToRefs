In development - a new workflow that takes fastq files and references as input and outputs the number of reads aligning to each reference.

Assumed directory structure:
.
├── Snakefile
├── data/
│   └── reads/
│       ├── sample1.fastq.gz
│       └── sample2.fastq.gz
├── references/
│   ├── refA.fasta
│   └── refB.fasta
├── alignments/
└── reports/