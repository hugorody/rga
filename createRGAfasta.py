#!/usr/bin/python3

from Bio import SeqIO

RGApredicted = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/RGAdb_beta2.tsv"
RGAorthologs = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/RGAdb_blast.tsv"
fasta_file = "/home/hugo/Dropbox/Esalq/references/AP85-441/monoploide_version/Sspon_monoploidy.protein.fasta"

RGAs = {}
with open(RGApredicted,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        if i[1] != "Non-RGA":
            RGAs[i[0]] = ''

with open(RGAorthologs,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        RGAs[i[0]] = ''

outRGAfasta = open("RGAS.fasta","w")
with open(fasta_file,"r") as set1:
    for i in SeqIO.parse(set1, "fasta"):
        if str(i.id) in RGAs:
            outRGAfasta.write(i.format("fasta"))
