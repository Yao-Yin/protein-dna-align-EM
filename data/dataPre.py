# test for get sequences

import requests, os, csv, json
import pandas as pd

df = pd.read_csv(r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\9606_pseudogenes.txt", sep="\t", header=0)
# dataPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\test_pg.txt"
# dataPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_450_pre_suf.txt"
# dataPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_pre_suf_10.txt"
# dataPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\blabla.txt"
dataPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10.txt"
errorPath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\test_pg_error.txt"
datafile = open(dataPath, "a+")
errorfile = open(errorPath, "a+")

valid_base_set = {"A", "T", "C", "G", 'a', 't', 'c', 'g'}
valid_aa_set = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'}
base_comp = {"A":"T", "C":"G", 'a':'t', 'c':'g', "T":"A", "G":"C", 't':'a', 'g':'c',}

def checkValid(dna_seq, pro_seq)->bool:
    for base in dna_seq:
        if base not in valid_base_set:
            return False
    for aa in pro_seq:
        if aa not in valid_aa_set:
            return False
    return True

def trans(dnaStr)->str:
    res = ''.join([base_comp[i] for i in dnaStr])
    return res[::-1]

def ucscJsonToFasta(originalText):
    jsonObject = json.loads(originalText)
    chr = jsonObject["chrom"]
    genome = jsonObject["genome"]
    start = jsonObject["start"]
    end = jsonObject["end"]
    dna = jsonObject["dna"]
    header = [genome, chr, start, end]
    res = ">" + "_".join([ str(x) for x in header])
    for i in range(len(dna)):
        if i % 60 == 0: res += "\n"
        res += dna[i].upper()
    return res


def getSequences(pseudogene_id, pseudogene_chr, pseudogene_start, pseudogene_end, strand, protein_id):
    ensemblServer = "https://rest.ensembl.org"
    proteinext = "/sequence/id/"+str(protein_id)+"?type=protein"
    ucscServer = "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;"
    dna_start_ext = "start="+str(pseudogene_start)+";"
    dna_end_ext = "end="+str(pseudogene_end)
    dna_chrom = "chrom=chr" + str(pseudogene_chr) + ";"
    try: 
        proteinSeqRequest = requests.get(ensemblServer+proteinext, headers={ "Content-Type" : "text/x-fasta"})
        print(proteinSeqRequest.text)
    except requests.HTTPError:
        print('Error with protein id code:', protein_id)
    try:
        dnaSeqRequest = requests.get(ucscServer + dna_chrom + dna_start_ext + dna_end_ext)
        print(ucscJsonToFasta(dnaSeqRequest.text))
    except requests.HTTPError:
        print("Error with protein id code;", pseudogene_id)

def getSmallTestData(pseudogene_id, pseudogene_chr, pseudogene_start, pseudogene_end, protein_id, strand, maxProLength = 200):
    ensemblServer = "https://rest.ensembl.org"
    proteinext = "/sequence/id/"+str(protein_id)+"?type=protein"
    ucscServer = "https://api.genome.ucsc.edu/getData/sequence?genome=hg19;"
    dna_start_ext = "start="+str(pseudogene_start)+";"
    dna_end_ext = "end="+str(pseudogene_end)
    dna_chrom = "chrom=chr" + str(pseudogene_chr) + ";"
    try: 
        proteinSeqRequest = requests.get(ensemblServer+proteinext, headers={ "Content-Type" : "text/x-fasta"})
        dnaSeqRequest = requests.get(ucscServer + dna_chrom + dna_start_ext + dna_end_ext)
    except requests.HTTPError:
        errorfile.write(str(pseudogene_id) + "\n")
        return
    curr_dna_json = json.loads(dnaSeqRequest.text)
    curr_protein_fasta = proteinSeqRequest.text
    curr_protein = curr_protein_fasta.split("\n")
    curr_protein_seq = "".join(curr_protein[1:])
    if len(curr_protein_seq) >= maxProLength :
        errorfile.write (str(pseudogene_id)+"; protein sequence is longer than "+ str(maxProLength) + "\n")
        return
    if not (proteinSeqRequest.ok and dnaSeqRequest.ok):
        errorfile.write (str(pseudogene_id)+" can not find protein.\n")
        return
    curr_pro_header = curr_protein[0]
    curr_dna_header = ">"+str(curr_dna_json["genome"]) + "_" + str(curr_dna_json["chrom"]) +"_"+ str(curr_dna_json["start"])+"_"+str(curr_dna_json["end"]) + "_" + strand
    curr_dna_seq = "".join(curr_dna_json["dna"].split("\n"))
    if checkValid(curr_dna_seq, curr_protein_seq):
        datafile.write(curr_pro_header+ "\n"+curr_protein_seq +"\n")
        datafile.write(curr_dna_header+"\n")
        if strand == '-':
            curr_dna_seq = trans(curr_dna_seq)
        datafile.write(curr_dna_seq + "\n")
    #datafile.write(">"+str(curr_dna_json["genome"]) + "_" + str(curr_dna_json["chrom"]) +"_"+ str(curr_dna_json["start"])+"_"+str(curr_dna_json["end"]) + "\n" + "".join(curr_dna_json["dna"].split("\n")) +"\n")
    

#PGOHUM00000250213	6	85967053	85967406	+	ENSP00000264258	0.95	0.75	3.0E-19	Processed	PF01198;
#PGOHUM00000242041	X	135929543	135930187	+	ENSP00000328551	1.0	0.791	0.0	Processed	PF00071;PF08477;
#getSequences("PGOHUM00000250213",6,85967053,85967406,"ENSP00000264258")

for idx, data in df.iterrows():
    if data["class"] == "Processed" and len(data["chr"]) <= 2 and int(data["end"]) - int(data["start"]) <= 500:
        getSmallTestData(data["id"], data["chr"], int(data["start"])-10, int(data["end"])+10, data["parent"], data["strand"], 300)


datafile.close()
errorfile.close()
