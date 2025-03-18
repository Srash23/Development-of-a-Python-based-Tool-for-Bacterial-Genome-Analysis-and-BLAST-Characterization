# Bacterial Genome ORF Analysis and BLAST Pipeline

## Introduction
This project develops a Python-based bioinformatics tool for the analysis of bacterial genomes, enabling:

1. Identification of Open Reading Frames (ORFs).
2. Translation of ORFs into protein sequences.
3. Molecular mass calculation of proteins.
4. BLAST (Basic Local Alignment Search Tool) analysis to find homologous sequences in public databases.

The sequence analyzed in this project is the complete synthetic genome of Mycoplasma mycoides JCVI-syn1.0, a milestone project by the Craig Venter Institute. Studying this genome provides insights into the minimal genetic requirements for life, contributing to biotechnology and synthetic biology.

## Why Is This Important?
Understanding bacterial genome structure is crucial for:

**🔬 Gene Identification –** Discovering essential genes in synthetic and natural bacterial genomes.

**🧬 Synthetic Biology & Biotechnology –** Engineering minimal genomes for research and industrial applications.

**🏥 Medicine –** Understanding pathogenic bacteria and developing targeted therapies.

**🌍 Evolutionary Biology –** Comparing synthetic genomes with natural bacteria to explore evolutionary relationships.

## Workflow
![Pipeline Workflow - visual selection](https://github.com/user-attachments/assets/f5a8d445-54b4-479c-913d-82633e5fe34c)

## Key Features
1. ORF Detection – Finds potential coding sequences in bacterial genomes.
2. Protein Translation – Converts ORFs into protein sequences.
3. Molecular Mass Calculation – Determines the size of translated proteins.
4. BLAST Analysis – Identifies homologous proteins in global databases.
5. Automated CSV Reports – Saves results for further study.

## Installation and Setup

**This project requires Python 3.7+ and the following dependencies:**

**Dependencies**

Install the required packages using:

pip install biopython requests

**Ensure you have:**

Genome sequence file (.fasta)

Internet access for BLAST queries.

## Key Insights

**1. Gene Identification –** Detected multiple ORFs, highlighting potential coding regions.

**2. Protein Characterization –** Translated sequences vary in molecular mass, useful for functional classification.

**3. Homology Detection –** BLAST analysis reveals evolutionary relationships with known bacterial proteins.

## Biological Significance

**1. Minimal Genome Analysis –** Identifies essential genes for bacterial survival.
   
**3. Comparative Genomics –** Compares synthetic bacterial genomes to natural species.
   
**5. Pathogen Research –** Can be extended to study virulent genes in pathogens.

## Conclusion

This Bacterial Genome ORF Analysis and BLAST Pipeline provides an automated bioinformatics workflow for identifying protein-coding genes, computing molecular properties, and searching for homologous sequences. The study of synthetic bacterial genomes like Mycoplasma mycoides JCVI-syn1.0 enhances our understanding of minimal life forms and accelerates biotechnological advancements.
