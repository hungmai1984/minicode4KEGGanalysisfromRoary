#!/usr/bin/env python3
import pandas as pd
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt

unique_files = {
    "VppirAB":"EggNogmapper_VppirAB.tabular",
    "VppirABTc":"EggNogmapper_VppirABTc.tabular",
    "VpTc":"EggNogmapper_VpTc.tabular"
}
background_file = "pan_genome_Eggnogmapper.tabular"
Kegg_names = {
    "map03070": "Bacterial secretion system",
    "map00051": "Fructose and mannose metabolism",
    "map00040": "Pentose and glucuronate interconversions",
    "map02030": "Bacterial chemotaxis",
    "map02020": "Two-component system",
    "map00360": "Phenylalanine metabolism",
    "map04112": "Cell cycle - Caulobacter",
    "map03030": "DNA replication",
    "map01120": "Microbial metabolism in diverse environments",
    "map02040": "Flagellar assembly",
    "map00520": "Amino sugar and nucleotide sugar metabolism"}
def kegg_pathway(target_file):
    df = pd.read_csv(target_file,sep='\t',comment='#',header=None)
    pathway_col =df[12].dropna().astype(str)
    counts={}
    total_kegg_genes=0
    for row in pathway_col:
        if row =='-' or not row: continue
        total_kegg_genes+=1
        pathways=[p.strip() for p in row.split(',') if p.strip().startswith('map')]
        for p in set(pathways):
            counts[p]=counts.get(p,0)+1
    return counts, total_kegg_genes
bg_counts,bg_total=kegg_pathway(background_file)
all_results=[]
for name,path in unique_files.items():
    st_counts, st_total = kegg_pathway(path)
    if st_total ==0: continue
    strain_kegg = []
    for kegg_id in sorted(bg_counts.keys()):
        k=st_counts.get(kegg_id,0)
        n=bg_counts[kegg_id]
        p_val = hypergeom.sf(k - 1, bg_total, n, st_total)
        strain_kegg.append({
            'Strain':name,
            'Pathway_ID':kegg_id,
            'Pathway_name':Kegg_names.get(kegg_id,kegg_id),
            'strain_count':k,
            'background_count':n,
            'Rich_Factor':k/n if n>0 else 0,
            'P_value':p_val})
    df_strain=pd.DataFrame(strain_kegg)
    _,q_val,_,_=multipletests(df_strain['P_value'], method='fdr_bh')
    df_strain['Q_value']=q_val
    all_results.append(df_strain)
final_df =pd.concat(all_results)
final_df.to_csv("keggg_enrichment_analysis.csv", index=False)
kegg_graph=final_df[final_df['P_value']<0.05].copy()
kegg_graph['-log10_q']=-np.log10(kegg_graph['Q_value'])
kegg_graph['Label']=kegg_graph['Pathway_name'] + "("+ kegg_graph['Strain'] + ")"
plt.figure(figsize=(10,8))
plt.scatter(kegg_graph['Rich_Factor'],kegg_graph['Label'],s=kegg_graph['strain_count']*14,c=kegg_graph['-log10_q'],cmap='RdYlBu_r', edgecolors='none', alpha=0.7)
plt.colorbar(label='-log10(Q_value)')
plt.xlabel('Rich_Factor')
plt.title('KEGG_Pathway_Enrichment_Analysis')
plt.tight_layout(rect=[0,0,0.85,1])
for size in [2,5,10]:
    plt.scatter([], [],s=size * 20,c='gray',edgecolors=None,label=f'{size} genes', alpha=0.6)
    plt.legend(title='Gene_Count',labelspacing=1.2,borderpad=1,frameon=True,loc='upper left',bbox_to_anchor=(1.3,1))
    plt.savefig('kegg_enrichment_plot.png')
print("Analysis complete-Bye bye")