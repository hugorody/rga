#!/usr/bin/python3
from Bio import SeqIO
import pandas as pd
################################################################################
# Inputs
coilsfile = "Sspon_monoploidy.protein.fasta.new.coils" #coils predictions output
tmhmmfile = "Sspon_monoploidy.protein.fasta.tmhmm" #tmhmm predictions output
interprofile = "Sspon_monoploidy.protein.fasta.tsv" #interpro predictions output
pfamfile = "Sspon_monoploidy.protein.fasta.pfamscan" #pfam_scan predctions output
phobiusfile = "Sspon_monoploidy.protein.fasta.phobius" #phobius predictions output

#outputs
#out_falsRGAs = open("RGAdb_false_beta2.tsv","w") #sequences not harboring RGA basic domains
out_trueRGAs = open("RGAdb_beta2_COORDINATES.tsv","w") #sequences harboring RGA basic domains
out_domains = open("RGAdb_domains_beta2_COORDINATES.txt","w")
################################################################################
# parse PfamScan predictions output
PfamDOMAINS = ["LRR", "NB-ARC", "TIR", "kinase", "Lzipper", "NB-LRR"]

# Rx_N: This entry represents the N-terminal domain found in many plant resistance
# proteins [1]. This domain has been predicted to be a coiled-coil (CC), however
# the structure shows that it adopts a four helical bundle fold

def standardize_namePFAM(name):
    if name != "NB-LRR" and "LRR" in name:
        name = "LRR"
    elif "kinase" in name:
        name = "kinase"
    elif "TIR" in name:
        name = "TIR"
    elif "Rx_N" in name:
        name = "CC"
    elif "Lzipper" in name:
        name = "LZ"
    else:
        name = name
    return name

pfamDBaccessions = {} #{Pfam accession:[domain1,domain2]}
pfamDBfeatures = {} #{seqID:[[pf_accession,pf_domain,start,end],[pf_accession,pf_domain,start,end]]}
with open(pfamfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "#" not in i:
            i = i.split()
            if i != []:
                seqID = i[0]
                start = i[1]
                end = i[2]
                pf_accession = i[5]
                pf_domain = standardize_namePFAM(i[6])

                if any(x in pf_domain for x in PfamDOMAINS):

                    #feed pfamDBfeatures
                    if seqID not in pfamDBfeatures:
                        pfamDBfeatures[seqID] = [[pf_accession,pf_domain,start,end]]
                    else:
                        addfeature = pfamDBfeatures[seqID]
                        checkpf = []
                        for u in addfeature:
                            checkpf.append(u[0])
                        if pf_accession not in checkpf:
                            addfeature.append([pf_accession,pf_domain,start,end])
                        pfamDBfeatures[seqID] = addfeature

                    #feed pfamDBaccessions
                    if pf_accession not in pfamDBaccessions:
                        pfamDBaccessions[pf_accession] = [pf_domain]
                    else:
                        addfeature = pfamDBaccessions[pf_accession]
                        if pf_domain not in addfeature:
                            addfeature.append(pf_domain)
                        pfamDBaccessions[pf_accession] = addfeature

print ("Sequences with domains predicted by PfamScan:",len(pfamDBfeatures))

################################################################################
# parse PHOBIUS predictions output
phobiusDBfeatures = {} #{sequence:"TMP"}
with open(phobiusfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "TM" not in i:
            i = i.split()
            sequence = i[0]
            if int(i[1]) > 1:
                start = i[3].split("-")[0][1:]
                end = i[3].split("-")[-1][:-1]
                phobiusDBfeatures[sequence] = [["Phobius","TM",start,end]]

print ("Sequences with domains predicted by PHOBIUS:",len(phobiusDBfeatures))

################################################################################
# parse COILS predictions output
def findCOIL(seq):
    coords_tmp = []
    coords = [] #list with start and end of coils
    index = 0
    for i in seq:
        if i == "x":
            passe = 1
        else:
            passe = 0
            if coords_tmp != []:
                coords.append([coords_tmp[0],coords_tmp[-1]])
                coords_tmp = []
        if passe == 1:
            coords_tmp.append(index)
        index += 1
    return coords

coilsDBfeatures = {}
with open(coilsfile,"r") as set1:
    for i in SeqIO.parse(set1, "fasta"):
        sequenceID = str(i.id) #all sequences in COILS output had CC predicted
        sequence = list(str(i.seq))
        if "x" in sequence:
            coils = findCOIL(sequence) # [[57, 77], [150, 170]]
            for u in coils:
                if sequenceID not in coilsDBfeatures:
                    coilsDBfeatures[sequenceID] = [["Coils","CC",u[0],u[1]]]
                else:
                    addcoil = coilsDBfeatures[sequenceID]
                    addcoil.append(["Coils","CC",u[0],u[1]])
                    coilsDBfeatures[sequenceID] = addcoil

print ("Sequences with domains predicted by COILS:",len(coilsDBfeatures))

################################################################################
# parse TMHMM2 predictions output
tmhmmDBfeatures_tmp = {}
with open(tmhmmfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "Number of predicted TMHs" in i:
            i = i.split()
            sequence = i[1]
            tmh = int(i[-1])
            if tmh >= 1: #more than one TM predicted by TMHMM2
                tmhmmDBfeatures_tmp[sequence] = ["TM"]

tmhmmDBfeatures = {}
with open(tmhmmfile,"r") as set1:
    for i in set1:
        i = i.rstrip()
        if "TMhelix" in i:
            i = i.split()
            sequenceID = i[0]
            start = i[3]
            end = i[4]
            if sequenceID in tmhmmDBfeatures_tmp:
                if sequenceID not in tmhmmDBfeatures:
                    tmhmmDBfeatures[sequenceID] = [["TMHMM2","TM",start,end]]
                else:
                    addval = tmhmmDBfeatures[sequenceID]
                    addval.append(["TMHMM2","TM",start,end])
                    tmhmmDBfeatures[sequenceID] = addval

print ("Sequences with domains predicted by TMHMM2:",len(tmhmmDBfeatures))

################################################################################
# parse INTERPRO predictions output

# key words expected to be found in INTERPRO output
interproFEATURES = ["Leucine-rich",
                    "LysM",
                    "Serine/Threonine",
                    "serine/threonine-protein",
                    "Serine-threonine/tyrosine-protein",
                    "kinase",
                    "NB-ARC",
                    "zipper",
                    "TIR"]

#dict used by standardize_nameINTERPRO func
featuresdic = {"Leucine-rich":"LRR",
                "LysM":"LysM",
                "Serine/Threonine":"STK",
                "serine/threonine-protein":"STK",
                "Serine-threonine/tyrosine-protein":"STTK",
                "kinase":"kinase",
                "NB-ARC":"NB-ARC",
                "zipper":"LZ",
                "TIR":"TIR"}

interproDBfeatures = {}
interproDBaccessions = {}
df1 = pd.read_csv(interprofile, sep='\t')
for i in range(len(df1)):
    seqID = df1.loc[i][0]
    database = df1.loc[i][3]
    accession = df1.loc[i][4]
    domain = df1.loc[i][5]
    start = df1.loc[i][6]
    end = df1.loc[i][7]
    ipr_accession = df1.loc[i][11]
    ipr_domain = df1.loc[i][12]
    if str(type(ipr_domain)) != "<class 'float'>" :
        my_intersect = set(ipr_domain.split(" ")).intersection(interproFEATURES)
        if my_intersect != set():
            for i in my_intersect:
                name_domain = featuresdic[i]
                if seqID not in interproDBfeatures:
                    interproDBfeatures[seqID] = [[ipr_accession,name_domain,start,end]]
                else:
                    addval = interproDBfeatures[seqID]
                    checkdomains = []
                    for u in addval: #create a list of domains
                        checkdomains.append(u[0])

                    if ipr_accession not in checkdomains:
                        addval.append([ipr_accession,name_domain,start,end])
                    interproDBfeatures[seqID] = addval

                #feed interproDBaccessions
                if ipr_accession not in interproDBaccessions:
                    interproDBaccessions[ipr_accession] = [name_domain]
                else:
                    addval = interproDBaccessions[ipr_accession]
                    if name_domain not in addval:
                        addval.append(name_domain)
                    interproDBaccessions[ipr_accession] = addval
'''
    else:
        if str(type(domain)) != "<class 'float'>" :
            my_intersect = set(domain.split(" ")).intersection(interproFEATURES)
            if my_intersect != set():
                print (database,accession,domain)
                for i in my_intersect:
                    name_domain = featuresdic[i]
                    if seqID not in interproDBfeatures:
                        interproDBfeatures[seqID] = [[accession,name_domain,start,end]]
                    else:
                        addval = interproDBfeatures[seqID]
                        checkdomains = []
                        for u in addval: #create a list of domains
                            checkdomains.append(u[0])

                        if accession not in checkdomains:
                            addval.append([accession,name_domain,start,end])
                        interproDBfeatures[seqID] = addval

                    #feed interproDBaccessions
                    if ipr_accession not in interproDBaccessions:
                        interproDBaccessions[ipr_accession] = [name_domain]
                    else:
                        addval = interproDBaccessions[ipr_accession]
                        if name_domain not in addval:
                            addval.append(name_domain)
                        interproDBaccessions[ipr_accession] = addval
'''
print ("Sequences with domains predicted by INTERPRO:",len(interproDBfeatures))

################################################################################
#DOMAINS DB

SequencesDomainsDB = {} # dict having all the domains predicted by all the algorithms
SequencesDomainsDBcoordinates = {}
def addDOMAINS(DICTIONARY):
    for i in DICTIONARY.items():
        if i[0] not in SequencesDomainsDB:
            my_domains = []
            for u in i[1]:
                my_domains.append(u[1])
            SequencesDomainsDB[i[0]] = my_domains

            my_coordinates = []
            for y in i[1]:
                my_coordinates.append(",".join(map(str,y)))
            SequencesDomainsDBcoordinates[i[0]] = my_coordinates
        else:
            addfeature = SequencesDomainsDB[i[0]]
            for u in i[1]:
                if u[1] not in addfeature:
                    addfeature.append(u[1])
            SequencesDomainsDB[i[0]] = addfeature

            addcoordinates = SequencesDomainsDBcoordinates[i[0]]
            for y in i[1]:
                addcoordinates.append(",".join(map(str,y)))
            SequencesDomainsDBcoordinates[i[0]] = addcoordinates



PHOBIUS = addDOMAINS(phobiusDBfeatures)
COILS = addDOMAINS(coilsDBfeatures)
TMHMM2 = addDOMAINS(tmhmmDBfeatures)
INTERPRO = addDOMAINS(interproDBfeatures)
PFAMSCAN = addDOMAINS(pfamDBfeatures)

################################################################################
# RGA structure definition

RGAbasicdomains = ["NB-ARC", "TIR", "LRR", "NB-LRR"]

RGAencodingfamily = {"TM-LRR":[["TM","LRR"],["TM","NB-LRR"],["TM","LysM"]],
                    "NBS-LRR":[["NB-LRR"],["NB-ARC","LRR"]]}

RGA_NBSLRR_subgroups = {"TNL":["TIR"],
                        "CNL":["CC"]}

RGAothersubgroups = {"TM-CC":["TM","CC"],
                    "TN":["TIR","NB-ARC"],
                    "CN":["CC","NB-ARC"],
                    "TL":["TIR","LRR"],
                    "TIR":["TIR"]}

newtrueRGAs = {}
for i in SequencesDomainsDB.items():
    sequence = i[0]
    sequence_domains = i[1]
    num_domains = len(set(i[1]))

    for j in RGAencodingfamily.items():
        encoding_family = j[0]
        domain_combinations = j[1] #list
        for u in domain_combinations:
            if len(set(u).intersection(sequence_domains)) == len(u):
                if encoding_family == "TM-LRR": #try to find a subgroup
                    RLK = ["kinase","STK","STTK"]
                    domain_interaction = set(RLK).intersection(sequence_domains)
                    if  domain_interaction != set(): #find TM-LRR subgroups through presence/ausence of kinase domains
                        subgroup_name = "RLK"
                        newtrueRGAs[sequence] = [encoding_family,subgroup_name,sequence_domains]
                    else:
                        subgroup_name = "RLP"
                        newtrueRGAs[sequence] = [encoding_family,subgroup_name,sequence_domains]

                elif encoding_family == "NBS-LRR":
                    for subgroup in RGA_NBSLRR_subgroups.items():
                        subgroup_name = subgroup[0]
                        subgroup_domains = subgroup[1]
                        for domain in subgroup_domains:
                            if domain in sequence_domains:
                                #encoding_family = "NBS-LRR"
                                newtrueRGAs[sequence] = [encoding_family,subgroup_name,sequence_domains]

                    if sequence not in newtrueRGAs:
                        newtrueRGAs[sequence] = [encoding_family,"Other",sequence_domains]

    if sequence not in newtrueRGAs:
        for l in RGAothersubgroups.items():
            subgroup_name = l[0]
            subgroup_domains = l[1] #list
            if len(subgroup_domains) == len(sequence_domains):
                if len(set(subgroup_domains).intersection(sequence_domains)) == len(subgroup_domains):
                    newtrueRGAs[sequence] = ["Other",subgroup_name,sequence_domains]
                    #print ("Other",subgroup_name,sequence,sequence_domains)

    #Last step verify basic RGA domains
    if sequence not in newtrueRGAs:
        if set(sequence_domains).intersection(RGAbasicdomains) != set():
            newtrueRGAs[sequence] = ["Other","Other",sequence_domains]
            #print (sequence,sequence_domains)

    if sequence not in newtrueRGAs:
        newtrueRGAs[sequence] = ["Non-RGA","NA",sequence_domains]

#write output
out_trueRGAs.write("Sequence ID\tRGA encoding family\tRGA Subgroup\tRGA Domains\tCoordinates (accession,domain,start,end)\n")
for i in newtrueRGAs.items():
    out_trueRGAs.write (i[0]+"\t"+i[1][0]+"\t"+i[1][1]+"\t"+", ".join(i[1][2])+"\t"+"; ".join(SequencesDomainsDBcoordinates[i[0]])+"\n")

################################################################################
# Print basic statistics

countCLASS = {}
for i in newtrueRGAs.items():
    clas = i[1][1]
    if clas not in countCLASS:
        countCLASS[clas] = 1
    else:
        addclas = countCLASS[clas]
        addclas += 1
        countCLASS[clas] = addclas

out_domains.write ("Number of "+ str(len(newtrueRGAs)) +" RGAs in classes:" + "\n")
for i in countCLASS.items():
    out_domains.write(i[0] + "\t" + str(i[1]) + "\n")

out_domains.write ("\nInterpro domains:\n")
for i in interproDBaccessions.items():
    out_domains.write(i[0] + "\t" + ",".join(map(str,i[1])) + "\n")

out_domains.write ("\nPfam domains:\n")
for i in pfamDBaccessions.items():
    out_domains.write(i[0] + "\t" + ",".join(i[1]) + "\n")
