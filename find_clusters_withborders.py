#!/usr/bin/python3

#Based on Cristian et al., 2016
#gene cluster
#- at least 3 NBS-LRR
#- having <9 other genes neighboring NBS-LRR genes
#- two neighboring NBS-LRR genes are <250kb apart

#gene supercluster :::::::: not being calculated
#- at least one NBS-LRR gene cluster, in addition to at least TWO other NBS-LRR genes
#- having <99 other genes neighboring NBS-LRR genes
#- two neighboring NBS-LRR genes are <2500kb apart

#1 kbp (kilo base pairs) = 1,000 bp (base pairs)
#1 Mbp (mega base pairs) = 1,000,000 bp
#1 Gbp (giga base pairs) = 1,000,000,000 bp

################################################################################
# INPTUS
from Bio import SeqIO

gff_file = "/home/hugo/Dropbox/Esalq/references/AP85-441/monoploide_version/Sspon_monoploidy.gff3"
RGApredicted = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/RGAdb_beta2_COORDINATES.tsv"
RGAorthologs = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/RGAdb_blast.tsv"

################################################################################
# CREATE GENE SORGHUM GENE COORDINATES DICTIONARY
chromosomes = ["1","2","3","4","5","6","7","8"]
chrcoordinates = {}
with open(gff_file,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "#" not in i:
            i = i.split("\t")
            if i[2] == "gene":
                chr = i[0].replace("Chr0","").replace("Chr","")
                start = int(i[3])
                end = int(i[4])
                sequenceID = i[-1].split(";")[0].replace("ID=","")
                if chr in chromosomes:
                    chrcoordinates[sequenceID] = [chr,start,end]

################################################################################
# CREATE RGA DICTIONARY
RGAS = {} # {"Specie RGA sequence ID": "RGA subgroup"}
with open(RGApredicted,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        sequenceID = i[0].split("|")[0].replace(".1","")
        if i[1] != "Non-RGA":
            RGAS[sequenceID] = i[2]

with open(RGAorthologs,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        sequenceID = i[0].split("|")[0].replace(".1","")
        RGAS[sequenceID] = i[6]
################################################################################
# FIND GENE CLUSTER

count_neighboring = 0
current_chr = []
tmp_cluster = []
gene_clusters = {}
count_cluster = 0

for i in sorted(chrcoordinates.items()): #starting reading gene by gene from 0 in chromosome 1
    sorghumID = i[0]
    chr = i[1][0]
    start = i[1][1]
    end = i[1][2]
    if current_chr == []:
        current_chr.append(chr)

    if chr == current_chr[-1]: #if sequence is in the same previous chromosome
        ################################################################### FIND
        if sorghumID in RGAS: #if sequence is RGA
            if tmp_cluster == []: #if tmp cluster is empty
                tmp_cluster.append(sorghumID) #start a new tmp cluster
                count_neighboring = 0
            else:
                if count_neighboring < 10:
                    atual_start = chrcoordinates[sorghumID][1]
                    previous_end = chrcoordinates[tmp_cluster[-1]][2]
                    if int(atual_start - previous_end) < 250000: #distance between two neighboring genes must be < 250kb
                        tmp_cluster.append(sorghumID)
                else:
                    if len(tmp_cluster) >= 3: #gene cluster must have at least 3 genes
                        gene_clusters[str(chr)+"_"+str(count_cluster)] = tmp_cluster #add cluster based on filters
                        count_cluster += 1
                        count_neighboring = 0 #start counting gene neighbors again
                        tmp_cluster = [] #start a new tmp cluster
                    else:
                        tmp_cluster = [] #start a new tmp cluster
                        tmp_cluster.append(sorghumID)
                        count_neighboring = 0
        else:
            count_neighboring += 1
        ################################################################### FIND

    else: #if new chromosome, start counting again
        count_neighboring = 0
        count_cluster = 0
        tmp_cluster = []

        ################################################################### FIND
        if sorghumID in RGAS: #if sequence is RGA
            if tmp_cluster == []: #if tmp cluster is empty
                tmp_cluster.append(sorghumID) #start a new tmp cluster
                count_neighboring = 0
            else:
                if count_neighboring < 10:
                    atual_start = chrcoordinates[sorghumID][1]
                    previous_end = chrcoordinates[tmp_cluster[-1]][2]
                    if int(atual_start - previous_end) < 250000: #distance between two neighboring genes must be < 250kb
                        tmp_cluster.append(sorghumID)
                else:
                    if len(tmp_cluster) >= 3: #gene cluster must have at least 3 genes
                        gene_clusters[str(chr)+"_"+str(count_cluster)] = tmp_cluster #add cluster based on filters
                        count_cluster += 1
                        count_neighboring = 0 #start counting gene neighbors again
                        tmp_cluster = [] #start a new tmp cluster
                    else:
                        tmp_cluster = [] #start a new tmp cluster
                        tmp_cluster.append(sorghumID)
                        count_neighboring = 0
        else:
            count_neighboring += 1
        ################################################################### FIND

    if chr not in current_chr:
        current_chr.append(chr)
################################################################################
# Statistics
selected_subgroups = ["RLK","RLP","CNL","TNL","TM-CC","CN"]
out_clusters = open("RGAclusters_withborders.tsv","w")
out_clusters.write("CHR\tCluster Number\tSpecies members\tCluster start\t"+
                    "Cluster end"+"\t"+"Cluster length"+"\t"+"Cluster size\tRGA subgroups (N)\t"+
                    "RGA number of subgroups\n")
for i in sorted(gene_clusters.items()):
    chr = i[0].split("_")[0]
    cluster_number = i[0].split("_")[1]
    start_cluster = chrcoordinates[i[1][0]][1]
    end_cluster = chrcoordinates[i[1][-1]][2]
    cluster_length = end_cluster - start_cluster
    startend = [start_cluster,end_cluster]
    sorghum_members = i[1]
    species_members = []
    RGAsubgroups = []
    for j in sorghum_members:
        species_members.append(j)

        if RGAS[j] in selected_subgroups:
            RGAsubgroups.append(RGAS[j])
        else:
            RGAsubgroups.append("Other")

    cluster_size = len(sorghum_members)

    RGAsubgroups_size = []
    has_subgroups = set(RGAsubgroups).intersection(selected_subgroups)

    if has_subgroups != set():
        for z in selected_subgroups:
            if z in RGAsubgroups:
                RGAsubgroups_size.append(z+"("+str(RGAsubgroups.count(z))+")")

    if "Other" in RGAsubgroups:
        RGAsubgroups_size.append("Other"+"("+str(RGAsubgroups.count("Other"))+")")

    out_clusters.write (chr+"\t"+cluster_number+
                        "\t"+",".join(species_members)+"\t"+"\t".join(map(str,startend))+
                        "\t"+str(cluster_length)+"\t"+str(cluster_size)+"\t"+
                        ",".join(RGAsubgroups_size)+"\t"+str(len(RGAsubgroups_size))+"\n")
