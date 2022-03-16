#! /usr/bin/env python


import glob
import pandas as pd
import numpy as np

print("packages")
## AA frequency
AA_usage_filename = glob.glob("*AA_content.csv")[0]
## Read depth for each contig
depth_file = glob.glob("*constats.txt")[0]

base_name = AA_usage_filename.replace('_AA_content.csv', "")


tax = pd.read_csv(glob.glob("*.out")[0], sep="\t", names = ["class", "gene_ID", "ncbiID", "length_score",
                                                           "taxonId", "accession", "match", "phylogeny"])

tax['Domain'] = tax.phylogeny.str.split(";", expand=True)[1]

tax['contig_ID'] = tax.gene_ID.str.replace('_[0-9]+$', "")

print("BP1")

tax_sub = tax[["gene_ID", "Domain", "contig_ID"]]
tax_sub.Domain = tax_sub.Domain.replace(np.nan, 'Unclassified', regex=True)

depth = pd.read_csv(depth_file, sep="\t")
depth = depth.rename(columns={"#ID":"contig_ID"})

keggs = pd.read_csv(glob.glob("*kofam_results.txt")[0], sep="\t", names = ["gene_ID", "KEGG"])
keggs['contig_ID'] = keggs.gene_ID.str.replace('_[0-9]+$', "")
keggs

AA_usage = pd.read_csv(AA_usage_filename)
AA_usage['length'] = AA_usage.End - AA_usage.Start + 1 # + 1 to count the start bp
AA_usage.longID = AA_usage.longID.str.replace('^>','')
AA_usage['contig_ID'] = AA_usage.longID.str.replace('_[0-9]+$', "")
AA_usage = AA_usage.drop(columns = ['Unnamed: 0', 'ID'])
AA_usage = AA_usage.rename(columns={"longID":"gene_ID"})


tax_depth = tax_sub.merge(depth[["contig_ID", "Avg_fold"]], on = 'contig_ID', how = 'left')
tax_depth_kegg = tax_depth.merge(keggs, on = ['gene_ID', 'contig_ID'], how = 'left')
tax_depth_kegg_aa = tax_depth_kegg.merge(AA_usage[['gene_ID', 'gc_cont', 'length', 'carbon',
                               'hydrogen', 'nitrogen', 'oxygen', 'sulfur', 'contig_ID']],
                     on = ['gene_ID', 'contig_ID'])

print("BP2")
tax_depth_kegg_aa = tax_depth_kegg_aa.assign(bp_total = tax_depth_kegg_aa["length"]* tax_depth_kegg_aa["Avg_fold"])
tax_depth_kegg_aa = tax_depth_kegg_aa.assign(gc_total = tax_depth_kegg_aa["bp_total"]* tax_depth_kegg_aa["gc_cont"])

tax_depth_kegg_aa[['carbon','hydrogen',
                       'nitrogen','oxygen','sulfur']] = tax_depth_kegg_aa[['carbon','hydrogen',
                       'nitrogen','oxygen','sulfur']].multiply(tax_depth_kegg_aa["Avg_fold"], axis="index")

taxa_output_file = base_name+"_genecounts_by_domain.csv"

taxa_counts_output = tax_depth_kegg_aa.groupby('Domain').agg(
    Domain_gene_counts=pd.NamedAgg(column="Avg_fold", aggfunc="sum"),
    GC = pd.NamedAgg(column="gc_total", aggfunc="sum"),
    BP = pd.NamedAgg(column="bp_total", aggfunc="sum"),
    AA_C = pd.NamedAgg(column="carbon", aggfunc="sum"),
    AA_H = pd.NamedAgg(column="hydrogen", aggfunc="sum"),
    AA_N = pd.NamedAgg(column="nitrogen", aggfunc="sum"),
    AA_O = pd.NamedAgg(column="oxygen", aggfunc="sum"),
    AA_S = pd.NamedAgg(column="sulfur", aggfunc="sum"),
    )

taxa_counts_output = taxa_counts_output.assign(GC_adj = taxa_counts_output['GC']/taxa_counts_output['BP'])
taxa_counts_output = taxa_counts_output.drop(columns = ['GC','BP'])

taxa_counts_output.to_csv(taxa_output_file)

kegg_taxa_out = tax_depth_kegg_aa.groupby(['Domain', 'KEGG']).agg(
    gene_counts=pd.NamedAgg(column="Avg_fold", aggfunc="sum"),
    GC = pd.NamedAgg(column="gc_total", aggfunc="sum"),
    BP = pd.NamedAgg(column="bp_total", aggfunc="sum"),
    )

kegg_taxa_out = kegg_taxa_out.assign(GC_adj = kegg_taxa_out['GC']/kegg_taxa_out['BP'])

kegg_taxa_out_long = kegg_taxa_out.pivot_table(index= "KEGG", columns="Domain",
                                               values=["gene_counts", "GC_adj"])

kegg_taxa_out_long.columns =[s1 + str(s2) for (s1,s2) in kegg_taxa_out_long.columns.tolist()]
kegg_taxa_out_long.reset_index(inplace=True)

kegg_output_filename = base_name+"_kegg_by_domain.csv"
kegg_taxa_out_long.to_csv(kegg_output_filename)




