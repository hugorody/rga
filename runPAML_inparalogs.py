#!/usr/bin/python3

from Bio import SeqIO
import subprocess
import glob,os
from Bio.Phylo.PAML import yn00

inparalogs_file = "/home/hugo/Dropbox/Esalq/references/RGA/OrthoMCL/pairs/inparalogs.txt"

R570_file = "/home/hugo/Dropbox/Esalq/references/MonoploidReference/single_tiling_path_cds.fna"
COMPGG_file = ""
AP85_file = "/home/hugo/Dropbox/Esalq/references/AP85-441/monoploide_version/Sspon_monoploidy.cds.fasta"
Sorghum_file = "/home/hugo/Dropbox/Esalq/references/SorghumPLAZA/cds.selected_transcript.sbi.fasta"

organism = "AP85-441" #Chose the reference: COMPGG, Sorghum, R570, or AP85-441

#catalog used by OrthoMCL
taxons = {"COMPGG":"Taxon0|",
            "Sorghum":"Taxon1|",
            "R570":"Taxon2|",
            "AP85-441":"Taxon3|"}

#Catalog of organisms and its FASTA files
files = {"COMPGG":COMPGG_file,
            "Sorghum":Sorghum_file,
            "R570":R570_file,
            "AP85-441":AP85_file}

#Dictionary containing all predited in-paralogs
inparalogs_dic = {}
with open(inparalogs_file,"r") as set1:
    for i in set1:
        i = i.rstrip().split("\t")
        inparalogs_dic[i[0]] = i[1]

#organism FASTA dictionary
organism_fasta = {}
with open(files[organism],"r") as set1:
    for i in SeqIO.parse(set1,"fasta"):
        if organism == "R570":
            seqID = str(i.id).replace("t","p")
        if organism == "AP85-441":
            seqID = str(i.id)
        if organism == "Sorghum":
            seqID = str(i.id)
        #print (seqID)
        organism_fasta[seqID] = str(i.seq)

#Create FASTA file for each pair of in-paralogs
#Align FASTA files of each pair of in-paralogs using the TranslatorX
for i in inparalogs_dic.keys():
    if taxons[organism] in i:
        copy1 = i.replace(taxons[organism],"")
        copy2 = inparalogs_dic[i].replace(taxons[organism],"")
        if copy1 in organism_fasta and copy2 in organism_fasta:
            seq1 = organism_fasta[copy1]
            seq2 = organism_fasta[copy2]

            namefasta = copy1+"_"+copy2+".fasta"
            namealnfasta = copy1+"_"+copy2+".aln.fasta"

            #create FASTA
            outFASTA = open(namefasta,"w")
            outFASTA.write(">"+copy1+"\n"+seq1[0:-3]+"\n>"+copy2+"\n"+seq2[0:-3]+"\n")
            outFASTA.close()
            #alingn FASTA
            #muscle = subprocess.Popen("muscle -in " + namefasta + " -out " + namealnfasta,shell=True)
            #muscle.wait()
            translatorx = subprocess.Popen("perl translatorX.pl -i "+namefasta+" -o "+namealnfasta+" -p M -t T -c 1",shell=True)
            translatorx.wait()

            os.remove(namefasta)

#remove unecessary data
#remove all html
html_infolder = glob.glob('*.html')
for x in html_infolder:
    os.remove(x)

#remove all log
log_infolder = glob.glob('*.log')
for x in log_infolder:
    os.remove(x)

#remove all aaseqs
aaseqs_infolder = glob.glob('*.aaseqs')
for x in aaseqs_infolder:
    os.remove(x)

#Convert FASTA to PHYLIP
for i in inparalogs_dic.keys():
    if taxons[organism] in i:
        copy1 = i.replace(taxons[organism],"")
        copy2 = inparalogs_dic[i].replace(taxons[organism],"")
        if copy1 in organism_fasta and copy2 in organism_fasta:
            seq1 = organism_fasta[copy1]
            seq2 = organism_fasta[copy2]

            namealnfasta = copy1+"_"+copy2+".aln.fasta"
            namealnphy = copy1+"_"+copy2+".aln.phy"
            nametranslator = copy1+"_"+copy2+".aln.fasta.nt_ali.fasta"

            #convert FASTA 2 PHYLIP
            fasta2phy = subprocess.Popen("perl Fasta2Phylip.pl " + nametranslator + " " + namealnphy,shell=True)
            fasta2phy.wait()

            #remove unecessary data
            os.remove(nametranslator)
            os.remove(copy1+"_"+copy2+".aln.fasta.aa_ali.fasta")
            os.remove(copy1+"_"+copy2+".aln.fasta.aaseqs.fasta")
            os.remove(copy1+"_"+copy2+".aln.fasta.nt12_ali.fasta")
            os.remove(copy1+"_"+copy2+".aln.fasta.nt1_ali.fasta")
            os.remove(copy1+"_"+copy2+".aln.fasta.nt2_ali.fasta")
            os.remove(copy1+"_"+copy2+".aln.fasta.nt3_ali.fasta")

#Run PAML yn00 for each PHYLIP alignment and create the final output table
output_yn00 = open(organism+"_yn00_finalresult.csv","w")
output_yn00.write("Seq1\tSeq2\tdS\tdS SE\tdN\tdN SE\n")
for i in inparalogs_dic.keys():
    if taxons[organism] in i:
        copy1 = i.replace(taxons[organism],"")
        copy2 = inparalogs_dic[i].replace(taxons[organism],"")
        if copy1 in organism_fasta and copy2 in organism_fasta:
            seq1 = organism_fasta[copy1]
            seq2 = organism_fasta[copy2]

            namealnphy = copy1+"_"+copy2+".aln.phy"
            print ("Runing yn00 for",namealnphy)
            #run PAML yn00 software
            yn = yn00.Yn00()
            yn.alignment = namealnphy
            yn.out_file = "yn_out.txt"
            yn.working_dir = "./"
            yn.commonf3x4 = 1
            yn.weighting = None
            yn.icode = None
            yn.ndata = None
            yn.verbose = None
            yn.run()
            result = yn00.read("yn_out.txt")[copy1]
            for uu in result.values():
                ds = str(uu["YN00"]["dS"])
                dsse = str(uu["YN00"]["dS SE"])
                dn = str(uu["YN00"]["dN"])
                dnse = str(uu["YN00"]["dN SE"])
                output_yn00.write (copy1+"\t"+copy2+"\t"+ds+"\t"+dsse+"\t"+dn+"\t"+dnse+"\n")

            os.remove(namealnphy)
