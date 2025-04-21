# ORF Discovery and Functional Annotation via BLAST

This repository presents a streamlined Python-based pipeline for identifying open reading frames (ORFs) from genomic DNA sequences, computing their molecular weights, and annotating them using NCBI’s remote BLAST API. The project is designed for quick, scriptable protein discovery and functional analysis from nucleotide data.

## Project Overview
The goal of this project is to:
- Parse and process genomic FASTA files
- Detect protein-coding ORFs on both strands (≥50 codons)
- Compute molecular weight of translated proteins
- Filter sequences for optimal BLAST performance
- Query the NCBI `nr` database for functional annotation
- Export structured results to a CSV

## Pipeline Summary

### 1. Input Handling
- Load FASTA files using Biopython’s `SeqIO`
- Handle single or multi-record nucleotide datasets

### 2. ORF Identification
- Search all 6 reading frames using standard codon table
- Filter sequences based on protein length

### 3. Molecular Mass Calculation
- Compute protein molecular weight (in kDa) using Biopython utilities

### 4. Functional Annotation
- Perform `blastp` remotely using `NCBIWWW.qblast`
- Extract top 5 hits and associated E-values

### 5. Output Generation
- Save ORF sequence, molecular weight, and BLAST annotations into a CSV

## Input/Output Description

| File            | Description                                      |
|-----------------|--------------------------------------------------|
| `sequence.fasta`| Input genomic or contig-level FASTA sequence     |
| `results.csv`   | Output file containing ORF, mass, and BLAST hits |

## Tech Stack

- **Language**: Python 3
- **Libraries**: Biopython, requests, csv
- **Tools**: NCBI BLAST API (`qblast`)

## License
MIT License
