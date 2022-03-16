#! /usr/bin/env python

import glob
import pandas as pd
import numpy as np

# Read in files
## Codon frequency
codon_input_file = glob.glob("*codon_frequency.csv")[0]
codon_frequency = pd.read_csv(codon_input_file)
## AA frequency
AA_usage = pd.read_csv(glob.glob("*AA_content.csv")[0])
## Read depth for each contig
depth = pd.read_csv(glob.glob("*constats.txt")[0], sep="\t")
## Contigs with bacteria 
Bact_contigs = pd.read_csv(glob.glob("*Bact_contigs.txt")[0], sep="\t", names=["contig_ID"])
## kegg numbers
keggs = pd.read_csv(glob.glob("*kofam_results.txt")[0], sep="\t", names = ["gene_ID", "KEGG"])

# Rename some columns, correct some labels 
depth = depth.rename(columns={"#ID":"contig_ID"})
AA_usage.longID = AA_usage.longID.str.replace('^>','')
AA_usage['contig_ID'] = AA_usage.longID.str.replace('_[0-9]+$', "")
AA_usage = AA_usage.drop(columns = ['Unnamed: 0', 'ID'])
AA_usage = AA_usage.rename(columns={"longID":"gene_ID"})
codon_frequency['contig_ID'] = codon_frequency.ID.str.replace('_[0-9]+$', "")
codon_frequency = codon_frequency.rename(columns={"ID":"gene_ID"})
Bact_contigs = Bact_contigs.drop_duplicates()
keggs['contig_ID'] = keggs.gene_ID.str.replace('_[0-9]+$', "")

# Merge the codon and AA frequency
joined_codon_AA = codon_frequency.merge(AA_usage, on=['gene_ID', 'contig_ID'], how = 'inner')

# Add in depth, remove non bacterial
joined_codon_AA_depth = joined_codon_AA.merge(depth, on = 'contig_ID', how = 'right')
joined_codon_AA_depth = Bact_contigs.merge(joined_codon_AA_depth, on = 'contig_ID', how = 'left')

# Add depth to kegg annotations
keggs = keggs.merge(depth, on = 'contig_ID', how = 'right')
# Isolate bacterial reads
bact_keggs = Bact_contigs.merge(keggs, on = 'contig_ID', how = 'left')

# Summarize kegg data by read depth

bact_kegg_sumarized = bact_keggs.groupby('KEGG').agg(
        kegg_sum=pd.NamedAgg(column="Avg_fold", aggfunc="sum")
    )


# Adjusting the counts of each codon and amino acid

joined_codon_AA_depth[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']] = joined_codon_AA_depth[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']].multiply(joined_codon_AA_depth["Avg_fold"], axis="index")

joined_codon_AA_depth[['A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*','carbon','hydrogen',
                       'nitrogen','oxygen','sulfur']] = joined_codon_AA_depth[['A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*','carbon','hydrogen',
                       'nitrogen','oxygen','sulfur']].multiply(joined_codon_AA_depth["Avg_fold"], axis="index")

# Sum up amino acids and codons
AA_Codon_sum = joined_codon_AA_depth[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG','A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*','carbon','hydrogen',
                       'nitrogen','oxygen','sulfur']].sum()

# Add the C:N mean for all reads
AA_Codon_sum = AA_Codon_sum.append(joined_codon_AA_depth[['c_to_n']].mean())
# Calculate and add the total C:N
AA_Codon_sum["Adjusted_c_to_n"] = AA_Codon_sum.carbon/AA_Codon_sum.nitrogen



# From the read depth, count the GCs and the bp's. A little check column in there too
depth['GC_total'] = depth['Length'] * depth['Ref_GC'] * depth['Avg_fold']
depth['bp_total'] = depth['Length'] * depth['Avg_fold']
depth['GC_check'] = depth['GC_total']/depth['bp_total']

# Sum those columns. This yields a depth adjusted GC content
depth_adjusted_GC_content = depth['GC_total'].sum()/depth['bp_total'].sum()

# filter out non-bacterial contigs 
bact_depth = Bact_contigs.merge(depth, on = 'contig_ID', how = 'left')
# Do the same calculation 
bact_depth_adjusted_GC_content = bact_depth['GC_total'].sum()/bact_depth['bp_total'].sum()

base_name = codon_input_file.replace('_codon_frequency.csv', "")


GC_output_file = open(base_name+"_GC_output.txt", "w")
GC_output_file.write("MGID\tGC\tGC_bact\tAA_CN\n")
GC_output_file.write(base_name+"\t"+str(depth_adjusted_GC_content)+"\t"+str(bact_depth_adjusted_GC_content)+
                     "\t"+str(AA_Codon_sum["Adjusted_c_to_n"]))
GC_output_file.close()

kegg_output_file = base_name+"_kegg_output.csv"

bact_kegg_sumarized.to_csv(kegg_output_file)


# get the total number of codons
total_codons = AA_Codon_sum[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']].sum()

# convert codon counts to frequencies
AA_Codon_sum[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']] = AA_Codon_sum[['TTT','TTC','TTA','TTG','CTT','CTC',
                       'CTA','CTG','ATT','ATC','ATA','ATG',
                       'GTT','GTC','GTA','GTG','TAT','TAC',
                       'TAA','TAG','CAT','CAC','CAA','CAG',
                       'AAT','AAC','AAA','AAG','GAT','GAC',
                       'GAA','GAG','TCT','TCC','TCA','TCG',
                       'CCT','CCC','CCA','CCG','ACT','ACC',
                       'ACA','ACG','GCT','GCC','GCA','GCG',
                       'TGT','TGC','TGA','TGG','CGT','CGC',
                       'CGA','CGG','AGT','AGC','AGA','AGG',
                       'GGT','GGC','GGA','GGG']].div(total_codons)

# Get total AA
total_AA = AA_Codon_sum[['A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*']].sum()
# Covert to frequencies 
AA_Codon_sum[['A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*']] = AA_Codon_sum[['A','C','D','E','F','G',
                       'H','I','K','L','M','N',
                       'P','Q','R','S','T','V',
                       'W','Y','*']].div(total_AA)

AA_codon_output_file = base_name+"_AA_codon_output.csv"
AA_Codon_sum.to_csv(AA_codon_output_file)

print("Made it to the end of the script with no errors")
