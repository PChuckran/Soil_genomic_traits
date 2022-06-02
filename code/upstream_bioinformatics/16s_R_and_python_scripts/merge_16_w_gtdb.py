#! /usr/bin/env python

import numpy as np
import re
import csv
import sys

input_file_name = sys.argv[1]
metadata = pd.read_csv("/scratch/pfc25/Marker_NEON/gtdb/bac120_metadata_r202_smm.tsv", delimiter="\t")
colNames = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
seq_info = pd.read_csv(input_file_name,  delimiter="\t", header=None)

seq_info.columns = colNames

seq_info["depth"] = seq_info["qseqid"].str.extract("(\d+)$").astype("int64")

seq_info["accession"] = seq_info["sseqid"].str.extract("^([^~]+)")

seq_info_merged = seq_info.merge(metadata, on=["accession"], how='left')

seq_info_merged["size_adjusted"] = seq_info_merged["genome_size"]*seq_info_merged["depth"]
seq_info_merged["gc_bases"] = seq_info_merged["size_adjusted"]*(seq_info_merged["gc_percentage"]/100)
seq_info_merged["coding_bases"] = seq_info_merged["size_adjusted"]*(seq_info_merged["coding_density"]/100)
sub_99 = seq_info_merged.loc[(seq_info_merged['pident'] >= 99)]
sub_95 = seq_info_merged.loc[(seq_info_merged['pident'] >= 95)]
sub_97 = seq_info_merged.loc[(seq_info_merged['pident'] >= 97)]

ags_99 = sub_99["size_adjusted"].sum()/sub_99["depth"].sum()
ags_97 = sub_97["size_adjusted"].sum()/sub_97["depth"].sum()
ags_95 = sub_95["size_adjusted"].sum()/sub_95["depth"].sum()
gc_99 = (sub_99["gc_bases"].sum()/sub_99["depth"].sum())/ags_99
gc_97 = (sub_97["gc_bases"].sum()/sub_97["depth"].sum())/ags_97
gc_95 = (sub_95["gc_bases"].sum()/sub_95["depth"].sum())/ags_95
coding_99 = (sub_99["coding_bases"].sum()/sub_99["depth"].sum())/ags_99
coding_97 = (sub_97["coding_bases"].sum()/sub_97["depth"].sum())/ags_97
coding_95 = (sub_95["coding_bases"].sum()/sub_95["depth"].sum())/ags_95

base = input_file_name.replace(".txt", "")
output_file_name = input_file_name.replace(".txt", "_summary.csv")
output_file = open("/scratch/pfc25/Marker_NEON/size_estimates/"+output_file_name, "w+")
wr = csv.writer(output_file, quoting=csv.QUOTE_ALL)

wr.writerow([base+"_99", ags_99, gc_99, coding_99, len(sub_99)])
wr.writerow([base+"_97", ags_97, gc_97, coding_97, len(sub_97)])
wr.writerow([base+"_95", ags_95, gc_95, coding_95, len(sub_95)])
output_file.close()
seq_info_merged.to_csv("/scratch/pfc25/Marker_NEON/full_output/"+base+"_full_out.csv")
