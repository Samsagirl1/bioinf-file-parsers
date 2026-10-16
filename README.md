# bioinformatics-file-parsers

A small package for working with bioinformatics files (FASTQ, FASTA, BLAST, GenBank). Made as part of a Python practicum in bioinformatics. All functions use only the Python standard library.

## Installation

Python 3.10 or newer is required. External libraries are not required.

```bash
git clone https://github.com/Samsagirl1/bioinf-file-parsers.git
cd bioinf-file-parsers
```

## Package structure

```
-/
|- README.md
|- main.py # utilities (imports, filter_fastq and other main functions)
|- helpers/
    |- fastq.py       # helper functions for filter_fastq
    |- bio_files.py   # helper functions for convert_multiline_fasta_to_oneline, parse_blast_output, select_genes_from_gbk_to_fasta
```

## Main functions

### `filter_fastq`
Filters reads/FASTQ sequences by GC content, nucleotide sequence length, and average read quality (Phred 33 scale). Filtering is performed for each read at the moment the input file is read. The result is saved to the "filtered" folder.

```python
from main import filter_fastq

filter_fastq(input_fastq = "example_fastq.fastq", output_fastq = "filtered_fastq.fastq", gc_bounds = (40, 60), length_bounds = (50, 150), quality_threshold = 20)
```

- `gc_bounds` and `length_bounds` accept an interval (min, max) or a single number, which will be considered the upper bound, with the lower bound defaulting to 0
- `quality_threshold` accepts a number - the minimum average read quality

### `convert_multiline_fasta_to_oneline` 
Takes a FASTA file as input and transforms each sequence so that each sequence (DNA/RNA/protein/...) fits on one line. As a result of the function, all transformed sequences will be saved to a new FASTA file.

```python
from main import convert_multiline_fasta_to_oneline

convert_multiline_fasta_to_oneline(input_fasta = "example_fasta.fasta", output_fasta = "output_fasta.fasta")
```

- The argument `output_fasta` is optional

### `parse_blast_output`
Takes from the BLAST algorithm report the first name from the Description column for each QUERY request and saves the unique names alphabetically, one per line.

```python
from main import parse_blast_output

parse_blast_output(input_file = "example_blast_result.txt", output_file = "parsed_blast_result.txt")
```

- The argument `output_file` is optional

### `select_genes_from_gbk_to_fasta`
Finds genes of interest in a GenBank file (.gbk) and saves to a FASTA file the amino acid sequences (the /translation qualifier) of their neighbors (flanking genes). The sequences of the genes of interest are not written to the output file.

```python
from main import select_genes_from_gbk_to_fasta

select_genes_from_gbk_to_fasta(input_gbk= "example_gbk.gbk", genes = ["sdhA", "sucC"], output_fasta = "output_fasta.fasta", n_before = 1, n_after = 1)
```
- `n_before`, `n_after` accept the number of genes before and after the gene of interest. By default, they are equal to 1
- `genes` accepts a gene of interest or a list of genes of interest
- `output_fasta` is optional

## Limitations

- The functions work only with files in the expected formats (see the example_data folder)
- `select_genes_from_gbk_to_fasta` reads only `CDS` features
- When neighboring genes are passed to the `select_genes_from_gbk_to_fasta` function, their sequences may be written to the output FASTA file as each other's neighbors

## Author
Valeria Kirenskaya
Telegram: @samsagirl
