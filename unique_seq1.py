#!/usr/bin/env python3
import pandas as pd
from Bio import SeqIO
import os
Input_files = {"presence_absence": "gene_presence_absence.Rtab","pan_genome": "pan_genome_reference.fa"}
def unique_sequences():
    for key, path in Input_files.items():
        if not os.path.exists(path):
            print(f"Error:{path} not found")
            return
    df = pd.read_csv(Input_files["presence_absence"], sep ='\t')
    all_strains = df.columns [1:]
    prefix=input("Enter prefix for srain (e.g., Vp_Tc): ").strip()
    group_cols = [col for col in all_strains if col.startswith(prefix)]
    other_cols = [col for col in all_strains if col not in group_cols]
    if not group_cols:
        print(f"no strains found")
        return
    print(f"{len(group_cols)} {prefix} vs {len(other_cols)} other_strains")
    group_gene = df[group_cols].sum(axis=1) == len(group_cols)
    other_gene = df[other_cols].sum(axis=1) == 0
    unique_genes= set(df[group_gene & other_gene]['Gene'])
    if not unique_genes:
        print("No genes found ")
        return
    output_file = f"unique_genes_{prefix}.fasta"
    print(f"Found {len(unique_genes)} genes. Extracting sequences ")
    counts = 0
    with open(output_file,"w") as final:
        for record in SeqIO.parse(Input_files["pan_genome"], "fasta"):
            matched_gene = None
            gene_name=record.id
            if gene_name in unique_genes:
                matched_gene = gene_name
            else:
                for gene in list(unique_genes):
                    if gene in record.description:
                        matched_gene = gene
                        break
            if matched_gene:
                SeqIO.write(record, final, "fasta")
                counts += 1
                unique_genes.discard(matched_gene)
        print(f"{prefix}: {counts} unique genes found and saved.successfully ")