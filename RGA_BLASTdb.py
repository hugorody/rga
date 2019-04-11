#!/usr/bin/python3

rgadb_file = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/RGAdb_beta2.tsv"
blast_besthit_file = "/home/hugo/Dropbox/Esalq/references/RGA/AP85-441/Monoploid/Sspon_monoploidy.protein.fasta_vs_128genes.blastp.besthit"
prgdb_vs_predictions = "/home/hugo/Dropbox/Esalq/references/RGA/PRGdb/PRGDB_vs_predictions_beta2.tsv"

outBLAST = open("RGAdb_blast.tsv","w")
outBLAST.write("Sequence ID\tPRGdb BestHit\tIdentity\tAlignment Length\tQuery coverage\tE-value\tRGA subgroup\tPredicted as\n")
################################################################################
RGAcandidates = {}
with open(rgadb_file,"r") as set1:
    for i in set1:
        i = i.rstrip()
        i = i.split("\t")
        gene = i[0]
        encodingfamily = i[1]
        RGAcandidates[gene] = encodingfamily

PRGDBdata = {}
with open(prgdb_vs_predictions,"r") as set3:
    for i in set3:
        i = i.rstrip()
        i = i.split("\t")
        PRGDBdata[i[0]] = [i[1],i[2]] #gene, rga subgroup

RGAsubgroups = {}
with open(blast_besthit_file,"r") as set2:
    for i in set2:
        i = i.rstrip()
        i = i.split("\t")
        identity = float(i[2])
        PRGDBid = "PRGDB" + i[1].split("_")[0]
        if i[0] not in RGAcandidates:
            if identity >= 40.0:
                outBLAST.write ("\t".join(i)+"\t"+PRGDBdata[PRGDBid][1]+"\t"+"BLASTp"+"\n")

                if PRGDBdata[PRGDBid][1] not in RGAsubgroups:
                    RGAsubgroups[PRGDBdata[PRGDBid][1]] = 1
                else:
                    addval = RGAsubgroups[PRGDBdata[PRGDBid][1]]
                    addval += 1
                    RGAsubgroups[PRGDBdata[PRGDBid][1]] = addval

        else: #i[0] is in RGA candidates
            #now verify if it was classified as RGA
            if RGAcandidates[i[0]] == "Non-RGA" and identity >= 40.0:
                outBLAST.write ("\t".join(i)+"\t"+PRGDBdata[PRGDBid][1]+"\t"+RGAcandidates[i[0]]+"\n")

                if PRGDBdata[PRGDBid][1] not in RGAsubgroups:
                    RGAsubgroups[PRGDBdata[PRGDBid][1]] = 1
                else:
                    addval = RGAsubgroups[PRGDBdata[PRGDBid][1]]
                    addval += 1
                    RGAsubgroups[PRGDBdata[PRGDBid][1]] = addval


for i in RGAsubgroups.items():
    print (i[0],i[1])
