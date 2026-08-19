#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import os
RTAB_FILE = 'gene_presence_absence.Rtab'
FASTA_FILE = 'pan_genome_reference.fa'
def extract_custom_sequences():
    if not os.path.exists(RTAB_FILE):
        print(f"Error: {RTAB_FILE} not found.")
        return
    df = pd.read_csv(RTAB_FILE, sep='\t')
    all_strains = df.columns[1:]
    print(f"Available strains in file: {len(all_strains)} total)")
    group_a_name = input("\nEnter prefix for Group A (e.g., Vp_Tc): ").strip()
    group_b_name = input("Enter prefix for Group B (e.g., Vp_pirAB): ").strip()
    exclude_name = input("Enter prefix to EXCLUDE (e.g., Vp-pirAB-Tc): ").strip()
    cols_a = [col for col in all_strains if group_a_name in col]
    cols_b = [col for col in all_strains if group_b_name in col]
    cols_exclude = [col for col in all_strains if exclude_name in col]
    if not cols_a or not cols_b:
        print("Error: no strains found!")
        return
    in_all_a = df[cols_a].all(axis=1)
    in_all_b = df[cols_b].all(axis=1)
    if cols_exclude:
        absent_in_exclude = df[cols_exclude].sum(axis=1) == 0
        mask = in_all_a & in_all_b & absent_in_exclude
    else:
        mask = in_all_a & in_all_b
    target_genes = set(df[mask]['Gene'])
    if not target_genes:
        print("\nNo genes found")
        return
    output_filename = f"Shared_{group_a_name}_{group_b_name}_no_{exclude_name}.fasta"
    print(f"Found {len(target_genes)} genes. Extracting from FASTA into {output_filename}...")
    extracted_count = 0
    with open(output_filename, "w") as output_file:
        for record in SeqIO.parse(FASTA_FILE, "fasta"):
            gene_name=record.id
            matched_gene = None          
            if gene_name in target_genes:
                matched_gene=gene_name
            else:
                for gene in list(target_genes):
                    if gene in record.description:
                        matched_gene=gene
                        break
            if matched_gene:                 
                SeqIO.write(record, output_file, "fasta")
                extracted_count += 1
                target_genes.discard(matched_gene)         
    print(f" {extracted_count} sequences written.successfully ")