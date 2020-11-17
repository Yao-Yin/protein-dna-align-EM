from math import exp

codon_table = {
    "A":{"GCT", "GCC", "GCA", "GCG",},
    "R":{"CGT", "CGC", "CGA", "CGG", "AGA", "AGG",},
    "N":{"AAT", "AAC",},
    "D":{"GAT", "GAC",},
    "C":{"TGT", "TGC",},
    "Q":{"CAA", "CAG",},
    "E":{"GAA", "GAG",},
    "G":{"GGT", "GGC", "GGA", "GGG",},
    "H":{"CAT", "CAC"},
    "I":{"ATT", "ATC", "ATA",},
    "L":{"CTT", "CTC", "CTA", "CTG", "TTA", "TTG",},
    "K":{"AAA", "AAG",},
    "M":{"ATG",},
    "F":{"TTT", "TTC",},
    "P":{"CCT", "CCC", "CCA", "CCG",},
    "S":{"TCT", "TCC", "TCA", "TCG", "AGT", "AGC",},
    "T":{"ACT", "ACC", "ACA", "ACG"},
    "W":{"TGG",},
    "Y":{"TAT", "TAC",},
    "V":{"GTT", "GTC", "GTA", "GTG"},
    "O":{"TAA", "TAG",},
    "U":{"TGA",},
    "*":{"TAA", "TAG", "TGA"}
}

codon_table_reverse = {}
for key, value in codon_table.items():
    for codon in value:
        codon_table_reverse[codon] = key

aa_list = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','O','U','*']

aa_order_dict = {}
for i in range(len(aa_list)):
    aa_order_dict[aa_list[i]] = i

triplet_list = ["TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", "TAC", "TAA", "TAG", "TGT", "TGC", "TGA", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG"]

outputFile = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\initProb\piProb.txt"

scoreMat, probMat = [[] for x in range(21)], [[] for x in range(21)]
for i in range(21):
    for j in range(64):
        curr_aa = aa_list[i]
        curr_tri = triplet_list[j]
        if curr_tri in codon_table[curr_aa]:
            scoreMat[i].append(1)
        else:
            scoreMat[i].append(-1)

tot = 0
for i in range(21):
    for j in range(64):
       tot += exp(scoreMat[i][j])

for i in range(21):
    for j in range(64):
        probMat[i].append(exp(scoreMat[i][j]) / tot)



with open(outputFile, "w") as f:
    for l in probMat:
        f.write(" ".join([str(x) for x in l]) + "\n")

