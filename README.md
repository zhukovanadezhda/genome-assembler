# De Bruijn Graph-Based Assembler

## Introduction

This repository contains an implementation of a de Bruijn graph-based assembler to assemble the genome of the enterovirus A71. This genome is particularly interesting due to its short length: 7408 nucleotides, linear, and non-segmented.

The provided fastq file was generated using the ART program [Huang 2011] with the following command:

```bash
art_illumina -i eva71.fna -ef -l 100 -f 20 -o eva71 -ir 0 -dr 0 -ir2 0 -dr2 0 -na -qL 41 -rs 1539952693
```

The reads have maximum quality (41) and contain no insertions. Only the reads corresponding to the 5' -> 3' strands are provided.

In the folder `debruijn-tp/data/`, you will find:
- `eva71.fna`: the genome of the virus of interest
- `eva71_plus_perfect.fq`: reads

## Dependency Installation

To set up the environment, run the following commands:

```bash
conda env create -f environment.yml
conda activate genome-assembler
```

## Usage

Clone the repository and navigate to the project folder:

```bash
git clone git@github.com:zhukovanadezhda/genome-assembler.git
cd genome-assembler
```

Run the assembler with:

```bash
python3 debruijn.py \
-i <input fastq file>  \ # single-end fastq file
-k <kmer size>         \ # optional, default is 21
-o <output file>         # file with the contigs
```

## Tests

You can test the program by running:

```bash
pytest --cov=debruijn
```

## Contact

If you have any questions, feel free to contact me via email: nadiajuckova@gmail.com
