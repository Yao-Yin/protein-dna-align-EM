#coding=utf-8
from numpy import ndarray, random


class Node:
    def __init__(self, name, seqIdx, distribution, single):
        self.name = name
        self.seqIdx = seqIdx
        self.distribution = distribution
        self.next = []
        self.prob = []
        self.single = single
    def transition(self, another, prob):
        self.next.append(another)
        self.prob.append(prob)
    def checkTransitions(self, eps=1e-5):
        tot = 0
        for p in self.prob:
            tot += p
        return len(self.next) == len(self.prob) and abs(1-tot) <= eps
    def checkDistribution(self, eps=1e-5):
        tot = 0
        if self.single:
            for p in self.distribution:
                tot += p
            return abs(tot-1) <= eps
        else:
            for p in self.distribution:
                for q in p:
                    tot += q
            return abs(tot-1) <= eps
    def updateTransition(self):
        pass
    def updateDistribution(self):
        pass

class PairHMM:
    def __init__(self):
        self.nodes = dict()
        self.seri = []
        self.iterTimes = 0
        self.convergedFlag = False
    def addnode(self, state):
        # state: an instance of Node
        self.nodes[state.name] = state
    def addtransition(self, state1, state2, prob):
        node1 = self.nodes[state1]
        node2 = self.nodes[state2]
    def unconverged(self):
        return self.iterTimes == 0 or not self.convergedFlag
            
    def forward(self):
        pass
    def backward(self):
        pass
    def BaumWelch(self, proSeq, dnaSeq, eps):
        pass
    def saveParameters(self, place):
        pass

class Tool:
    def __init__(self):
        self.codonTable = {
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
        }
        self.baseOrd = {
            'T':0,
            'C':1,
            'A':2,
            'G':3,
        }
        self.Ordbase = ['T', 'C', 'A', 'G']
        self.RcodonTable = dict()
        for aa, codons in self.codonTable:
            for codon in codons:
                self.RcodonTable[codon] = aa
        
    def encodeAA(self, aa):
        return ord(aa) - ord('A')

    def decodeAA(self, idx):
        return chr(ord('A') + idx)

    def encodeDNATriplet(self, dna):
        res = 0
        for base in dna:
            res *= 4
            res += self.baseOrd[base]
        return res

    def decodeDNATriplet(self, code):
        third = self.Ordbase[code % 4]
        code //= 4
        second = self.Ordbase[code % 4]
        code //= 4
        first = self.Ordbase[code % 4]
        return first + second + third


