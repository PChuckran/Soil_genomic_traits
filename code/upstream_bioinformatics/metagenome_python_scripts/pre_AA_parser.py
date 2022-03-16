#! /usr/bin/env python


import sys
import csv
import pandas as pd


input_file = sys.argv[1]
output_file_name = input_file.strip(".faa")+"_AA_content.csv"
testline_raw = str()
c = str()

output_file = open(output_file_name, "w")
wr = csv.writer(output_file, quoting=csv.QUOTE_ALL)


i =0
# Create dictionary which will count AA
aa_counts = {"A":0, "C":0, "D":0, "E":0, "F":0, "G":0,
             "H":0, "I":0, "K":0, "L":0, "M":0, "N":0,
            "P":0, "Q":0, "R":0, "S":0, "T":0, "V":0, "W":0, "Y":0, "*":0}
colnames = ["longID", "Start", "End", "findout",
                "ID", "partial", "start_type", "rbs_motif",
                "rbs_spacer", "gc_cont"]
# Put together all column names
colnames.extend(list(aa_counts.keys()))
wr.writerow(colnames)
with open(input_file, 'r') as list_paths:
    # For each line, identify if it is a label or string of AA and parse accordingly
    for line in list_paths:
        print(i)
        if line.startswith(">"):
            # If it's the very first line, go straight to parse
            if i == 0:
                next
            else:
                # Marks the end of a sequence and the start of a new one
                # Append the last sequence to the DF
                testline.extend(list(aa_counts.values()))
                wr.writerow(testline)
                # Re-initialize dictionary so everything is back to zero
                aa_counts = {"A":0, "C":0, "D":0, "E":0, "F":0, "G":0,
                             "H":0, "I":0, "K":0, "L":0, "M":0, "N":0,
                            "P":0, "Q":0, "R":0, "S":0, "T":0, "V":0, "W":0, "Y":0, "*":0}
            i += 1
            testline = line.strip("\n").replace(" # ", ";").split(';')
        else:
            # Parse AA sequence 
            seq = line.strip("\n")
            for aa in seq:
                if aa in aa_counts:
                    aa_counts[aa] = aa_counts[aa] + 1
                else:
                    aa_counts[aa] = 1 
    # Append the last line to DF
    testline.extend(list(aa_counts.values()))
    wr.writerow(testline)
output_file.close()

data = pd.read_csv(output_file_name)

# Calculate total number of atoms for Carbon, Hydrogen
# Nitrogen, Oxygen, and Sulfur, for each amino acid

data["carbon"] = ((data["A"]*3) + (data["C"]*3)+ (data["D"]*4)+
                  (data["E"]*5) + (data["F"]*9)+ (data["G"]*2)+
                  (data["H"]*6) + (data["I"]*6)+ (data["K"]*6)+
                  (data["L"]*6) + (data["M"]*5)+ (data["N"]*4)+
                  (data["P"]*5) + (data["Q"]*5)+ (data["R"]*6)+
                  (data["S"]*3) + (data["T"]*4)+ (data["V"]*5)+
                  (data["W"]*11) + (data["Y"]*9)
                 )

data["hydrogen"] = ((data["A"]*7) + (data["C"]*7)+ (data["D"]*7)+
                  (data["E"]*9) + (data["F"]*11)+ (data["G"]*5)+
                  (data["H"]*9) + (data["I"]*13)+ (data["K"]*14)+
                  (data["L"]*13) + (data["M"]*11)+ (data["N"]*8)+
                  (data["P"]*9) + (data["Q"]*10)+ (data["R"]*14)+
                  (data["S"]*7) + (data["T"]*9)+ (data["V"]*11)+
                  (data["W"]*12) + (data["Y"]*11)
                 )

data["nitrogen"] = ((data["A"]*1) + (data["C"]*1)+ (data["D"]*1)+
                  (data["E"]*1) + (data["F"]*1)+ (data["G"]*1)+
                  (data["H"]*3) + (data["I"]*1)+ (data["K"]*2)+
                  (data["L"]*1) + (data["M"]*1)+ (data["N"]*2)+
                  (data["P"]*1) + (data["Q"]*2)+ (data["R"]*4)+
                  (data["S"]*1) + (data["T"]*1)+ (data["V"]*1)+
                  (data["W"]*2) + (data["Y"]*1)
                 )

data["oxygen"] = ((data["A"]*2) + (data["C"]*2)+ (data["D"]*4)+
                  (data["E"]*4) + (data["F"]*2)+ (data["G"]*2)+
                  (data["H"]*2) + (data["I"]*2)+ (data["K"]*2)+
                  (data["L"]*2) + (data["M"]*2)+ (data["N"]*3)+
                  (data["P"]*2) + (data["Q"]*3)+ (data["R"]*2)+
                  (data["S"]*3) + (data["T"]*3)+ (data["V"]*2)+
                  (data["W"]*2) + (data["Y"]*3)
                 )

data["sulfur"] = ((data["A"]*0) + (data["C"]*1)+ (data["D"]*0)+
                  (data["E"]*0) + (data["F"]*0)+ (data["G"]*0)+
                  (data["H"]*0) + (data["I"]*0)+ (data["K"]*0)+
                  (data["L"]*0) + (data["M"]*1)+ (data["N"]*0)+
                  (data["P"]*0) + (data["Q"]*0)+ (data["R"]*0)+
                  (data["S"]*0) + (data["T"]*0)+ (data["V"]*0)+
                  (data["W"]*0) + (data["Y"]*0)
                 )

# Calculate the C to N ratio of each gene call
data["c_to_n"] = data["carbon"]/data["nitrogen"]
data["c_to_n"] = pd.to_numeric(data["c_to_n"])

# Strip some weird formatting on the gc content
data["gc_cont"] = data["gc_cont"].str.strip("gc_cont=")
data["gc_cont"] = pd.to_numeric(data["gc_cont"])

data.to_csv(output_file_name)

