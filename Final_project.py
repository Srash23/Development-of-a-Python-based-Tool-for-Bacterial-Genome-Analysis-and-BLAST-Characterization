#!/usr/bin/env python
# coding: utf-8

# In[55]:


pip install biopython


# In[56]:


pip install biopython requests


# In[66]:


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data import CodonTable
from Bio.SeqUtils import molecular_weight
from Bio.Blast import NCBIWWW, NCBIXML
# Only import Requests if you're not using Bio.Blast for BLAST queries
import requests
import csv


# In[67]:


def read_genome(filename):
    """Reads a genome from a FASTA file."""
    with open(filename, 'r') as file:
        records = list(SeqIO.parse(file, "fasta"))
    print(f"Loaded {len(records)} records from {filename}.")
    print("__________________________________________________________________________________________________________")
    return records


# In[68]:


def find_orfs(sequence, min_protein_length=50):
    """Finds ORFs in a given sequence."""
    table = CodonTable.unambiguous_dna_by_id[11]  # Using the standard bacterial, archaeal and plant plastid code
    orfs = []
    print("------------------------------------")
    print("[1] Finding ORFs:")
    print("------------------------------------")
    for strand, nuc in [(+1, sequence), (-1, sequence.reverse_complement())]:
        for frame in range(3):
            length = 3 * ((len(nuc)-frame) // 3)  # Multiple of three
            for pro in nuc[frame:frame+length].translate(table).split("*"):
                if len(pro) >= min_protein_length:
                    orfs.append(pro)
                    print(f"Found ORF: Length {len(pro)} codons, Sequence {str(pro)[:100]}...")  # Print first 100 bases
    print("------------------------------------")
    print(f"Total ORFs found: {len(orfs)} , printing first 100 codons")
    print("__________________________________________________________________________________________________________")
    return orfs


# In[69]:


def calculate_molecular_mass(proteins):
    """Calculates the molecular mass of protein sequences."""
    print("------------------------------------")
    print("[2] Calculated molecular masses for all proteins in kDa.")
    print("------------------------------------")
    masses = []
    for protein in proteins:
        mass = molecular_weight(protein, seq_type='protein') / 1000.0  # Convert to kDa
        masses.append(mass)
        print(f"Protein sequence: {protein[:10]}... (length: {len(protein)}) | Molecular mass: {mass:.2f} kDa")
    print("__________________________________________________________________________________________________________")
    return masses


# In[70]:


def filter_sequences_for_blast(proteins, max_length=1000):
    """Filters proteins for BLAST analysis based on length."""
    filtered = [p for p in proteins if len(p) < max_length]
    print("[3] Filtered proteins for BLAST analysis (keeping it <1000 for ease of the analysis)")
    print("------------------------------------")
    print(f"Filtered down to {len(filtered)} proteins for BLAST.")
    print("------------------------------------")
    for i, protein in enumerate(filtered, 1):
        print(f"Protein {i}: {protein[:10]}... (length: {len(protein)})")  # Shows the first 10 amino acids and length    
    print("__________________________________________________________________________________________________________")
    return filtered


# In[71]:


def perform_blast(protein_sequences):
    """Performs BLAST search for a list of protein sequences."""
    print("[4] Performing BLAST")
    print("------------------------------------")
    print("Selected top 5 hits")
    blast_results = []
    for sequence in protein_sequences:
        result_handle = NCBIWWW.qblast("blastp", "nr", sequence)
        blast_record = NCBIXML.read(result_handle)
        hits = [(alignment.title, hsp.expect) for alignment in blast_record.alignments for hsp in alignment.hsps]
        blast_results.append(hits[:5])  # Top 5 hits
    print("Performed BLAST analysis.")
    print("__________________________________________________________________________________________________________")
    return blast_results


# In[72]:


def save_results_to_csv(results, filename="results.csv"):
    """Saves results to a CSV file."""
    print("[5] Saving data to a CSV file")
    print("------------------------------------")
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(["ORF", "MolecularMass (kDa)", "BLAST Hits"])
        for result in results:
            writer.writerow(result)
    print(f"Results saved to {filename}.")
    print("__________________________________________________________________________________________________________")


# In[73]:


# Ask for the genome file path
genome_file = "sequence.fasta"
genome_records = read_genome(genome_file)


# In[74]:


orfs = find_orfs(genome_records[0].seq)
proteins = [str(orf) for orf in orfs]  # Assuming ORFs are Seq objects
masses = calculate_molecular_mass(proteins)
filtered_proteins = filter_sequences_for_blast(proteins)
blast_results = perform_blast(filtered_proteins[:5])  # Limiting to 5 

compiled_results = list(zip(orfs, masses, blast_results))

save_results_to_csv(compiled_results)


# In[ ]:


# Used ChatGPT for debugging

