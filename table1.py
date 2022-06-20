#!/usr/bin/python3

from Bio import SeqIO
import re

RGApredicted = "/home/hugo/Dropbox/Esalq/Sunflower/Ha412-Eugene/RGAs/RGAdb_2020.tsv"
RGAorthologs = "/home/hugo/Dropbox/Esalq/Sunflower/Ha412-Eugene/RGAs/RGAdb_blast.tsv"

selected_subgroups = ["RLK","RLP","CNL","TNL","TM-CC","CN","TN","TIR"]
RGAS = {} # {sequence ID: RGAsubgroup}

with open(RGApredicted,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        if i[0] != 'Sequence ID':
            #sequenceID = i[0].replace("-RA","").replace("-RB","")
            sequenceID = i[0]
            if i[1] != "Non-RGA":
                if i[2] in selected_subgroups:
                    subgroup = i[2]
                else:
                    subgroup = "Other"
                RGAS[sequenceID] = subgroup

with open(RGAorthologs,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        if i[0] != 'Sequence ID':
            sequenceID = i[0]
            #sequenceID = re.sub("-R[A-Z]$","",i[0])
            if i[6] in selected_subgroups:
                subgroup = i[6]
            else:
                subgroup = "Other"
            RGAS[sequenceID] = subgroup

subgroups = {}
for i in RGAS.items():
    if i[1] not in subgroups:
        subgroups[i[1]] = 1
    else:
        subgroups[i[1]] = subgroups[i[1]] + 1

for i in subgroups.items():
    print (i)

print ("Total:",len(RGAS))
