import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from math import log
from pandas import *

class DataTool:
    def __init__(self):
        self.baseList = ['T', 'C', 'A', 'G']
        self.baseMap = {}
        for i in range(len(self.baseList)):
            self.baseMap[self.baseList[i]] = i
        
        self.aaList = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V', '*']
        self.aaMap = {}
        for i in range(len(self.aaList)):
            self.aaMap[self.aaList[i]] = i
        
        self.tripletList = [self.decodeTriplet(x) for x in range(64)]
        self.tripletMap = {}
        for i in range(len(self.tripletList)):
            self.tripletMap[self.tripletList[i]] = i
        
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
            "O":{"TAA", "TAG",},
            "U":{"TGA",},
            "*":{"TAA", "TAG", "TGA"}
        }

    def encodeAA(self, aaChr):
        return self.aaMap[aaChr]
    
    def decodeAA(self, aaCode):
        return self.aaList[aaCode]

    def encodeBase(self, baseChr):
        return self.baseMap[baseChr]

    def decodeBase(self, baseCode):
        return self.baseList[base]

    def encodeTriplet(self, TripletStr):
        return self.tripletMap[TripletStr]

    def decodeTriplet(self, TripletCode):
        first = (TripletCode >> 4) & 3 
        second = (TripletCode >> 2) & 3
        third = (TripletCode) & 3  
        return self.baseList[first] + self.baseList[second] + self.baseList[third]

    def checkCodon(self, aaCode, codenCode):
        pass


def plot_table(row, col, vals):
    """
    函数功能: 绘制二维表格，草图如下:
        -----------------
            |col1 |col2 |
        -----------------
        row1|value|value|
        -----------------
        row2|value|value|
        -----------------
    输入：
        row:string,(N)           #['row1', 'row2']
        col:string,(M)           #['col1', 'col2']
        vals:np, (N,M)          
    """
    
    R, C = len(row), len(col)
    idx = Index(row)
    df = DataFrame(np.random.randn(R, C), index=idx, columns=col)
    
    # 根据行数列数设置表格大小
    figC, figR = 2.25*C, R
    fig = plt.figure(figsize=(figC, figR))
    
    # 设置fig并去掉边框
    ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    the_table=plt.table(cellText=vals, rowLabels=df.index, colLabels=df.columns, colWidths = [0.01]*vals.shape[1], rowLoc='center', loc='center',cellLoc='center')
    the_table.set_fontsize(15)
    
    # 伸缩表格大小常数
    the_table.scale(figR/R*2 ,figC/C*1.5)
    plt.savefig("testdoc.jpg")

PiInputFilePath = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\visualization\testdoc"
dt = DataTool()

with open(PiInputFilePath, 'r') as f:
    # 21 * 64 matrix
    PiMatrix = []
    #for i in range(20):
    PiMatrix.append([log(float(x)) for x in f.readline().rstrip().split("\t")])
    #print(PiMatrix)
    row = dt.aaList
    col = dt.tripletList
    PiMatrixArray = np.array(PiMatrix)
    PiMatrixArray.resize((21, 64))
    print(PiMatrixArray)
    #print(PiMatrixArray)
    plot_table(np.array(row), np.array(col), PiMatrixArray)

    