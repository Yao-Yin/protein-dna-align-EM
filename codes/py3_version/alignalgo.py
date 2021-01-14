import numpy as np
from math import log, exp
import random as rd
class DataTool:
    def __init__(self):
        self.baseList = ['T', 'C', 'A', 'G']
        self.baseMap = {}
        for i in range(len(self.baseList)):
            self.baseMap[self.baseList[i]] = i
        self.aadistribution = [0.074, 0.042, 0.044, 0.059, 0.033, 0.058, 0.037, 0.074, 0.029, 0.038, 0.076, 0.072, 0.018, 0.04, 0.05, 0.081, 0.062, 0.013, 0.033, 0.068]
        self.basedistribution = [0.303, 0.199, 0.303, 0.195]
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
            "*":{"TAA", "TAG", "TGA"}
        }

    def encodeAA(self, aaChr):
        return self.aaMap[aaChr]
    
    def decodeAA(self, aaCode):
        return self.aaList[aaCode]

    def encodeBase(self, baseChr):
        return self.baseMap[baseChr]

    def decodeBase(self, baseCode):
        return self.baseList[baseCode]

    def encodeTriplet(self, TripletStr):
        return self.tripletMap[TripletStr]

    def decodeTriplet(self, TripletCode):
        first = (TripletCode >> 4) & 3 
        second = (TripletCode >> 2) & 3
        third = (TripletCode) & 3  
        return self.baseList[first] + self.baseList[second] + self.baseList[third]
    def genRandDNA(self, seq_len):
        d_code = [x for x in range(4)]
        dna = [ self.genRndnum(d_code, self.basedistribution) for x in range(seq_len)]
        return "".join([self.decodeBase(x) for x in dna])
    def genRndnum(self, num_list, pro_list):
        x = rd.uniform(0, 1)
        # 累积概率
        cum_pro = 0.0
        # 将可迭代对象打包成元组列表
        for number, number_pro in zip(num_list, pro_list):
            cum_pro += number_pro
            if x < cum_pro:
                # 返回值
                return number
    def genRandPro(self, seq_len):
        p_code = [x for x in range(20)]
        pro = [ self.genRndnum(p_code, self.aadistribution) for x in range(seq_len)]
        return "".join([self.decodeAA(x) for x in pro])

class Tracer():
    def __init__(self, prevx, prevy, prevMat, value):
        self.x = prevx
        self.y = prevy
        self.prevMat = prevMat
        self.nx = -1
        self.ny = -1
        self.val = value
        self.nextMat = None
    def setnext(self, nx, ny, nextMat):
        self.nx = nx
        self.ny = ny
        self.nextMat = nextMat
    
    def setprev(self, prevx, prevy, prevMat):
        self.x = prevx
        self.y = prevy
        self.prevMat = prevMat
    
    def setValue(self, v):
        self.val = v

ws = [[0 for i in range(1000)] for j in range(1000)]
wp = [[0 for i in range(1000)] for j in range(1000)]
#dt = DataTool()
#print(sum(dt.aadistribution))
#exit(0)

class Parameters:
    def __init__(self):
        # Transition probs:
        # omega_i, omega_d, gamma, alpha_i, alpha_d,
        # delta_i, beta_i, epsilon_i,
        # delta_d, beta_d, epsilon_d,
        # Emission probs:
        # phi: 21
        # psi: 4
        # pi: 21 * 64
        # OverallProbs
        # different Scores:
        # SxY, a_i, b_i, a_d, b_d, f_i, f_d, g_i, g_d,
        self.omega_i = 0.0
        self.omega_d = 0.0
        self.gamma = 0.0
        self.alpha_i = 0.0
        self.alpha_d = 0.0
        self.delta_i = 0.0
        self.beta_i = 0.0
        self.epsilon_i = 0.0
        self.delta_d = 0.0
        self.beta_d = 0.0
        self.epsilon_d = 0.0
        self.phi = [0.0 for x in range(21)]
        self.psi = [0.0 for x in range(4)]
        self.pi = [[0.0 for x in range(64)] for y in range(21)]
        self.prob = 0.0
        self.s = [[0.0 for x in range(64)] for y in range(21)]
        self.a_i = 0.0
        self.b_i = 0.0
        self.a_d = 0.0
        self.b_d = 0.0
        self.f_i = 0.0
        self.g_i = 0.0
        self.f_d = 0.0
        self.g_d = 0.0
        self.dt = DataTool()

    def setTransitions(self, firstline):
        tran_probs = [float(x) for x in firstline.rstrip().split('\t')]
        self.omega_i = tran_probs[0]
        self.omega_d = tran_probs[1]
        self.gamma = tran_probs[2]
        self.alpha_i = tran_probs[3]
        self.alpha_d = tran_probs[4]
        self.delta_i = tran_probs[6]
        self.beta_i = tran_probs[7]
        self.epsilon_i = tran_probs[8]
        self.delta_d = tran_probs[9]
        self.beta_d = tran_probs[10]
        self.epsilon_d = tran_probs[11]

    def setPhi(self, phiStr):
        self.phi = [float(x) for x in phiStr.rstrip().split('\t')]

    def setPsi(self, psiStr):
        self.psi = [float(x) for x in psiStr.rstrip().split('\t')]

    def setPi(self, piStr):
        curr_pi = [float(x) for x in piStr.rstrip().split('\t')]
        for i in range(21):
            for j in range(64):
                self.pi[i][j] = curr_pi[i*64+j]

    def reNormalize(self):
        tot_phi, tot_psi, tot_pi = 0.0, 0.0, 0.0
        for i in range(len(self.phi)):
            tot_phi += self.phi[i]
        for i in range(len(self.psi)):
            tot_psi += self.psi[i]
        for i in range(len(self.pi)):
            for j in range(len(self.pi[0])):
                tot_pi += self.pi[i][j]
        for i in range(len(self.phi)):
            self.phi[i] /= tot_phi
        for i in range(len(self.psi)):
            self.psi[i] /= tot_psi
        for i in range(21):
            for j in range(64):
                self.pi[i][j] /= tot_pi
    def display(self):
        print("a_i:", self.a_i)
        print("b_i:", self.b_i)
        print("a_d:", self.a_d)
        print("b_d:", self.b_d)
        print("f_i:", self.f_i)
        print("f_d:", self.f_d)
        print("g_i:", self.g_i)
        print("g_d:", self.g_d)
    def get_new_score(self, x, y):
        return self.pi[x][y]/(sum([i for i in self.pi[x]])*sum([self.pi[i][y] for i in range(21)]))

    def fixTest(self):
        for x in range(21):
            for y in range(64):
                self.s[x][y] = log(self.get_new_score(x, y)*self.gamma/(self.omega_d*self.omega_i**3))

    def getScore(self):
        # different Scores:
        # SxY, a_i, b_i, a_d, b_d, f_i, f_d, g_i, g_d,
        self.a_i = log(self.alpha_i*(1-self.beta_i)/self.beta_i)
        self.a_d = log(self.alpha_d*(1-self.beta_d)/self.beta_d)
        self.b_i = (1/3)*log((self.beta_i*self.delta_i*self.epsilon_i)/self.omega_i)
        self.b_d = (1/3)*log((self.beta_d*self.delta_d*self.epsilon_d)/self.omega_d)
        self.f_i = log(((1-self.delta_i)/(1-self.beta_i))) + (1/3)*log((self.beta_i**2)/(self.beta_i*self.epsilon_i))
        self.g_i = log((1-self.delta_i)/(1-self.beta_i)) + (1/3)*log(self.beta_i*self.delta_i/(self.epsilon_i**2))
        self.f_d = log((1-self.delta_d)/(1-self.beta_d)) + (1/3)*log((self.beta_d**2)/(self.delta_d*self.epsilon_d*(self.omega_d**2))) - 2*log(self.omega_i)
        self.g_d = log((1-self.epsilon_d)/(1-self.beta_d)) + (1/3)*log((self.beta_d*self.delta_d)/(self.omega_d*(self.epsilon_d**2))) - log(self.omega_i)
        for x in range(21):
            for y in range(64):
                self.s[x][y] = log(self.getSubstitionScore(x, y)*self.gamma/(self.omega_d*self.omega_i**3))
    def getNewScore(self, stat):
        for x in range(21):
            for y in range(64):
                self.s[x][y] = log(self.getNewSub(x, y, stat)*self.gamma/(self.omega_d*self.omega_i**3))
    def getNewSub(self, x, y, stat):
        return self.pi[x][y]/(self.phi[x]*stat[y])
    def getSubstitionScore(self, x, y):
        # return SxY
        return self.pi[x][y]/(self.phi[x]*self.psi[y & 3]*self.psi[(y >> 2) & 3]*self.psi[(y >> 4) & 3])
    def getS(self, x, Y):
        return 
    def alignScore(self, pro_seq, dna_seq):
        n = len(pro_seq)
        m = len(dna_seq)
        pro_num = [-1 if i == 0 else self.dt.encodeAA(pro_seq[i-1]) for i in range(n + 1)]
        dna_num = [-1 if i == 0 else self.dt.encodeBase(dna_seq[i-1]) for i in range(m + 1)]
        dna_tri = [self.dt.encodeTriplet(dna_seq[i-1:i+2]) if i > 0 and i + 2 <= m else -1 for i in range(m+1)]
        #print(pro_num)
        #print(dna_num)
        #print(dna_tri)  
        X = np.full((n+1, m+1), -1000.0, dtype=np.double)
        Y = np.full((n+1, m+1), -1000.0, dtype=np.double)
        Z = np.full((n+1, m+1), -1000.0, dtype=np.double)
        X[0][0] = 0
        for i in range(0,n+1):
            for j in range(0, m+1):
                y1 = (Y[i-1][j-2] + self.b_d + self.f_d) if i > 0 and j > 1 else -1000.0
                y2 = (Y[i-1][j-1] + 2*self.b_d + self.g_d) if i > 0 and j > 0 else -1000.0
                y3 = (Y[i-1][j] + 3*self.b_d) if i > 0 else -1000.0
                z1 = (Z[i][j-1] + self.b_i + self.f_i) if j > 0 else -1000.0
                z2 = (Z[i][j-2] + 2*self.b_i + self.g_i) if i > 0 and j > 1 else -1000.0
                z3 = (Z[i][j-3] + 3*self.b_i) if j > 2 else -1000.0
                w = max(0, y1, y2, y3, z1, z2, z3)
                w = max(w, X[i-1][j-3]) if i > 0 and j > 2 else w
                
                #print("s", i, j , w)
                if i > 0 and j > 0:
                    X[i][j] = w + self.s[pro_num[i]][dna_tri[j]] 
                elif i == 0 and j == 0:
                    X[i][j] = 0
                else:
                    X[i][j] = -1000.0
                if j > 0:  
                    Y[i][j] = max(w+self.a_d, y3)
                if i > 0:
                    Z[i][j] = max(w + self.a_i, z3) 
                # if i < 2 and j < 2: print(i, j, Y[i][j])
                # ws[i][j] = Y[i][j]
        res = -10000
        for i in range(n+1):
            for j in range(m+1):
                res = max(res, X[i][j])
        return res
    def alignScore2(self, pro_seq, dna_seq):
        n = len(pro_seq)
        m = len(dna_seq)
        pro_num = [-1 if i == 0 else self.dt.encodeAA(pro_seq[i-1]) for i in range(n + 1)]
        dna_num = [-1 if i == 0 else self.dt.encodeBase(dna_seq[i-1]) for i in range(m + 1)]
        dna_tri = [self.dt.encodeTriplet(dna_seq[i-1:i+2]) if i > 0 and i + 2 <= m else -1 for i in range(m+1)]
        #print(pro_num)
        #print(dna_num)
        #print(dna_tri)  
        X = np.full((n+1, m+1), -1000.0, dtype=np.double)
        Y = np.full((n+1, m+1), -1000.0, dtype=np.double)
        w = np.full((n+1, m+1), -1000.0, dtype=np.double)
        z1 = np.full((n+1, m+1), -1000.0, dtype=np.double)
        z2 = np.full((n+1, m+1), -1000.0, dtype=np.double)
        z3 = np.full((n+1, m+1), -1000.0, dtype=np.double)
        X[0][0] = 0
        for i in range(0,n+1):
            for j in range(0, m+1):
                y1 = (Y[i-1][j-2] + self.b_d + self.f_d) if i > 0 and j > 1 else -1000.0
                y2 = (Y[i-1][j-1] + 2*self.b_d + self.g_d) if i > 0 and j > 0 else -1000.0
                y3 = (Y[i-1][j] + 3*self.b_d) if i > 0 else -1000.0
                
                z1[i][j] = max(z1[i][j-3] + 3*self.b_i if j > 3 else -1000.0, w[i][j-1] + self.b_i + self.f_i + self. a_i) if i > 0 and j > 1 else -1000.0
                z2[i][j] = max(z2[i][j-3] + 3*self.b_i if j > 3 else -1000.0, w[i][j-2] + 2*self.b_i + self.g_i + self. a_i) if i > 0 and j > 2 else -1000.0
                z3[i][j] = max(z3[i][j-3] + 3*self.b_i if j > 3 else -1000.0, w[i][j-3] + 3*self.b_i + self. a_i) if i > 0 and j > 3 else -1000.0
                w[i][j] = max(0, y1, y2, y3)
                w[i][j] = max(w[i][j], X[i-1][j-3]) if i > 0 and j > 2 else w[i][j]
                k = max(w[i][j], z1[i][j], z2[i][j], z3[i][j])
                #print("s", i, j , w)
                if i > 0 and j > 0:
                    X[i][j] = k + self.s[pro_num[i]][dna_tri[j]] 
                elif i == 0 and j == 0:
                    X[i][j] = 0
                else:
                    X[i][j] = -1000.0
                if i > 0:  
                    Y[i][j] = max(k+self.a_d, y3)
                ws[i][j] = z3[i][j]
                # if i > 0:
                #    Z[i][j] = max(w + self.a_i, z3) 
                # if i < 2 and j < 2: print(i, j, Y[i][j])
                # ws[i][j] = Y[i][j]
        res = -10000
        for i in range(n+1):
            for j in range(m+1):
                res = max(res, X[i][j])
        #for i in range(4):
        #    for j in range(4):
        #        print(i, j, y3[i][j])
        # #print("######")
        return res
    
    def alignPattern(self, pro_seq, dna_seq):
        n = len(pro_seq)
        m = len(dna_seq)
        pro_num = [-1 if i == 0 else self.dt.encodeAA(pro_seq[i-1]) for i in range(n + 1)]
        dna_num = [-1 if i == 0 else self.dt.encodeBase(dna_seq[i-1]) for i in range(m + 1)]
        dna_tri = [self.dt.encodeTriplet(dna_seq[i-1:i+2]) if i > 0 and i + 2 <= m else -1 for i in range(m+1)]
        #print(pro_num)
        #print(dna_num)
        #print(dna_tri)  
        X = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        Y = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        Z = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        w = [[Tracer(-1, -1, None, 0.0) for x in range(m+1)] for y in range(n+1)]
        y1 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        y2 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        y3 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z1 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z2 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z3 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        X[0][0].val = 0
        for i in range(0, n+1):
            for j in range(0, m+1):
                if i > 0 and j > 2:
                    y1[i][j] = Tracer(i-1, j-2, "Y", Y[i-1][j-2].val + self.b_d + self.f_d) 
                if i > 0 and j > 1:
                    y2[i][j] = Tracer(i-1, j-1, "Y", Y[i-1][j-1].val + 2*self.b_d + self.g_d) 
                if i > 0:
                    y3[i][j] = Tracer(i-1, j , "Y", Y[i-1][j].val + 3*self.b_d) 
                if j > 0:
                    z1[i][j] = Tracer(i, j-1, "Z", Z[i][j-1].val + self.b_i + self.f_i) 
                if i > 0 and j > 1:
                    z2[i][j] = Tracer(i, j-2, "Z", Z[i][j-2].val + 2*self.b_i + self.g_i) 
                if j > 2:
                    z3[i][j] = Tracer(i, j-3, "Z", Z[i][j-3].val + 3*self.b_i)  
                w_val = max(0, y1[i][j].val, y2[i][j].val, y3[i][j].val, z1[i][j].val, z2[i][j].val, z3[i][j].val, (X[i-1][j-3].val if i > 0 and j > 2 else 0))
                
                w[i][j].setValue(w_val)
                if w_val == 0:
                    pass
                elif i > 0 and j > 2 and w_val == X[i-1][j-3].val:
                    w[i][j].setprev(i-1, j-3, "X")
                elif w_val == y3[i][j].val:
                    w[i][j].setprev(i, j, "y3")
                elif w_val == z3[i][j].val:
                    w[i][j].setprev(i, j, "z3")
                elif w_val == y2[i][j].val:
                    w[i][j].setprev(i, j, "y2")
                elif w_val == z2[i][j].val:
                    w[i][j].setprev(i, j, "z2")
                elif w_val == y1[i][j].val:
                    w[i][j].setprev(i, j, "y1")
                elif w_val == z1[i][j].val:
                    w[i][j].setprev(i, j, "z1")
                X[i][j].setprev(i, j, "w")
                if i > 0 and j > 0:
                    X[i][j].setValue(w[i][j].val + self.s[pro_num[i]][dna_tri[j]]) 
                elif i == 0 and j == 0:
                    X[i][j].setValue(0)
                if i > 0:
                    if w[i][j].val + self.a_d >= y3[i][j].val:
                        Y[i][j] = Tracer(i, j, "w", w[i][j].val + self.a_d)
                    else:
                        Y[i][j] = Tracer(i, j, "y3", y3[i][j].val)
                if j > 0:
                    if w[i][j].val + self.a_i >= z3[i][j].val:
                        Z[i][j] = Tracer(i, j, "w", w[i][j].val + self.a_i)
                    else:
                        Z[i][j] = Tracer(i, j, "z3", z3[i][j].val)
                # wp[i][j] = Y[i][j].val
        res = -10000
        # Backtracing
        curr_i, curr_j = -1, -1
        for i in range(n+1):
            for j in range(m+1):
                if X[i][j].val > res:
                    res = X[i][j].val
                    curr_i, curr_j  = i, j
        # Backtracing
        # part1 : handle the end
        head_pro, tail_pro = "", ""
        head_dna, tail_dna = "", ""
        align_dna = []
        align_pro = []
        tail_pro = pro_seq[curr_i:]
        tail_dna = dna_seq[curr_j+2:]
        #part 2: get align region
        # case 1: amino acid X align to DNA triplet yyy
        # #X#
        # yyy
        # case 2: Deletion
        # &X&
        # --y
        # case 3: Insertion
        # --
        # yy
        curr_tracer = X[curr_i][curr_j]
        curr_type = "X"
        while curr_tracer.prevMat != None:
            if curr_type == "X":
                # R_i align Q_j
                align_pro.append("#" + self.dt.decodeAA(pro_num[curr_i]) + "#")
                align_dna.append(self.dt.decodeTriplet(dna_tri[curr_j]))
            if curr_type == "y1":
                # 3n - 2 deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = dna_seq[curr_j-3:curr_j-1]
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("-" + curr_ins_dna)
            if curr_type == "y2":
                # 3n - 1 deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = dna_seq[curr_j-2:curr_j-1]
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("--" + curr_ins_dna)
            if curr_type == "y3":
                # 3n deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = ""
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("---")
            if curr_type == "z1":
                # 3n - 2 insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                curr_ins_dna = dna_seq[curr_j-2]
                curr_del_pro = "-"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)            
            if curr_type == "z2":
                # 3n - 1 insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                curr_ins_dna = dna_seq[curr_j-3:curr_j-1]
                curr_del_pro = "--"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)
            if curr_type == "z3":
                # 3n insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                curr_ins_dna = dna_seq[curr_j-4:curr_j-1]
                curr_del_pro = "---"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)
            nx_type = curr_tracer.prevMat
            
            if nx_type == None:
                break
            curr_type = nx_type
            curr_i = curr_tracer.x
            curr_j = curr_tracer.y
            #print(nx_type, curr_i, curr_j)
            if nx_type == "w":
                curr_tracer = w[curr_i][curr_j]
            if nx_type == "X":
                curr_tracer = X[curr_i][curr_j]
            if nx_type == "y1":
                curr_tracer = y1[curr_i][curr_j]
            if nx_type == "y2":
                curr_tracer = y2[curr_i][curr_j]
            if nx_type == "y3":
                curr_tracer = y3[curr_i][curr_j]
            if nx_type == "z1":
                curr_tracer = z1[curr_i][curr_j]
            if nx_type == "z2":
                curr_tracer = z2[curr_i][curr_j]
            if nx_type == "z3":
                curr_tracer = z3[curr_i][curr_j]
            if nx_type == "Y":
                curr_tracer = Y[curr_i][curr_j]
            if nx_type == "Z":
                curr_tracer = Z[curr_i][curr_j]
        head_pro = pro_seq[:curr_i-1]
        head_dna = dna_seq[:curr_j-1]
        align_region_pro = "".join(x for x in align_pro[::-1])
        align_region_dna = "".join(x for x in align_dna[::-1])
        print(align_region_pro)
        print(align_region_dna)
        return head_pro , head_dna , align_region_pro, align_region_dna, tail_pro, tail_dna, res

    def alignPattern2(self, pro_seq, dna_seq):
        n = len(pro_seq)
        m = len(dna_seq)
        pro_num = [-1 if i == 0 else self.dt.encodeAA(pro_seq[i-1]) for i in range(n + 1)]
        dna_num = [-1 if i == 0 else self.dt.encodeBase(dna_seq[i-1]) for i in range(m + 1)]
        dna_tri = [self.dt.encodeTriplet(dna_seq[i-1:i+2]) if i > 0 and i + 2 <= m else -1 for i in range(m+1)]
        #print(pro_num)
        #print(dna_num)
        #print(dna_tri)  
        X = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        Y = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        Z = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        w = [[Tracer(-1, -1, None, 0.0) for x in range(m+1)] for y in range(n+1)]
        k = [[Tracer(-1, -1, None, 0.0) for x in range(m+1)] for y in range(n+1)]
        y1 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        y2 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        y3 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z1 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z2 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        z3 = [[Tracer(-1, -1, None, -1000.0) for x in range(m+1)] for y in range(n+1)]
        X[0][0].val = 0
        for i in range(0, n+1):
            for j in range(0, m+1):
                if i > 0 and j > 1:
                    y1[i][j] = Tracer(i-1, j-2, "Y", Y[i-1][j-2].val + self.b_d + self.f_d) 
                if i > 0 and j > 0:
                    y2[i][j] = Tracer(i-1, j-1, "Y", Y[i-1][j-1].val + 2*self.b_d + self.g_d) 
                if i > 0 :
                    y3[i][j] = Tracer(i-1, j , "Y", Y[i-1][j].val + 3*self.b_d) 

                if i > 0 and j > 1:
                    if j > 3 and z1[i][j-3].val + 3*self.b_i >= w[i][j-1].val + self.f_i + self.b_i + self.a_i:
                        z1[i][j] = Tracer(i, j-3, "z1", z1[i][j-3].val + 3*self.b_i)
                    else:
                        z1[i][j] = Tracer(i, j-1, "w", w[i][j-1].val + self.b_i + self.f_i + self.a_i) 

                if i > 0 and j > 2:
                    if j > 3 and  z2[i][j-3].val + 3*self.b_i >= w[i][j-2].val + self.g_i + 2*self.b_i + self.a_i:
                        z2[i][j] = Tracer(i, j-3, "z2", z2[i][j-3].val + 3*self.b_i)
                    else:
                        z2[i][j] = Tracer(i, j-2, "w", w[i][j-2].val + 2*self.b_i + self.g_i + self.a_i) 

                if i > 0 and j > 3:
                    if z3[i][j-3].val + 3*self.b_i >= w[i][j-3].val + self.a_i + 3*self.b_i:
                        z3[i][j] = Tracer(i, j-3, "z3", z3[i][j-3].val + 3*self.b_i)
                    else:
                        z3[i][j] = Tracer(i, j-3, "w", w[i][j-3].val + 3*self.b_i + self.a_i) 
                
                w_val = max(0, y1[i][j].val, y2[i][j].val, y3[i][j].val, (X[i-1][j-3].val if i > 0 and j > 2 else 0))
                w[i][j].setValue(w_val)
                if w_val == 0:
                    pass
                elif w_val == y3[i][j].val:
                    w[i][j] = Tracer(i, j, "y3", w_val)
                elif w_val == y2[i][j].val:
                    w[i][j] = Tracer(i, j, "y2", w_val)
                elif w_val == y1[i][j].val:
                    w[i][j] = Tracer(i, j, "y1", w_val)
                elif i > 0 and j > 2 and w_val == X[i-1][j-3].val:
                    w[i][j] = Tracer(i-1, j-3, "X", w_val)
                k_val = max(w_val, z1[i][j].val, z2[i][j].val, z3[i][j].val)
                k[i][j].setValue(k_val)
                
                if k_val == 0:
                    pass
                elif i > 0 and j > 2 and k_val == X[i-1][j-3].val:
                    k[i][j].setprev(i-1, j-3, "X")
                elif k_val == y3[i][j].val:
                    k[i][j].setprev(i, j, "y3")
                elif k_val == z3[i][j].val:
                    k[i][j].setprev(i, j, "z3")
                elif k_val == y2[i][j].val:
                    k[i][j].setprev(i, j, "y2")
                elif k_val == z2[i][j].val:
                    k[i][j].setprev(i, j, "z2")
                elif k_val == y1[i][j].val:
                    k[i][j].setprev(i, j, "y1")
                elif k_val == z1[i][j].val:
                    k[i][j].setprev(i, j, "z1")
                if i > 0 and j > 0:
                    X[i][j].setprev(i, j, "k")
                    X[i][j].setValue(k_val + self.s[pro_num[i]][dna_tri[j]]) 
                elif i == 0 and j == 0:
                    X[i][j].setValue(0)
                if i > 0:
                    if k[i][j].val + self.a_d >= y3[i][j].val:
                        Y[i][j] = Tracer(i, j, "k", k[i][j].val + self.a_d)
                    else:
                        Y[i][j] = Tracer(i, j, "y3", y3[i][j].val)
                wp[i][j] = z3[i][j].val
        res = -10000
        curr_i, curr_j = -1, -1
        for i in range(n+1):
            for j in range(m+1):
                if X[i][j].val > res:
                    res = X[i][j].val
                    curr_i, curr_j  = i, j
        # Backtracing
        # part1 : handle the end
        head_pro, tail_pro = "", ""
        head_dna, tail_dna = "", ""
        align_dna = []
        align_pro = []
        tail_pro = pro_seq[curr_i:]
        tail_dna = dna_seq[curr_j+2:]
        #part 2: get align region
        # case 1: amino acid X align to DNA triplet yyy
        # #X#
        # yyy
        # case 2: Deletion
        # &X&
        # --y
        # case 3: Insertion
        # --
        # yy
        curr_tracer = X[curr_i][curr_j]
        curr_type = "X"
        while curr_tracer.prevMat != None:

            if curr_type == "X":
                # R_i align Q_j
                align_pro.append("#" + self.dt.decodeAA(pro_num[curr_i]) + "#")
                align_dna.append(self.dt.decodeTriplet(dna_tri[curr_j]))
            if curr_type == "y1":
                # 3n - 2 deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = dna_seq[curr_j-3:curr_j-1]
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("-" + curr_ins_dna)
            if curr_type == "y2":
                # 3n - 1 deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = dna_seq[curr_j-2:curr_j-1]
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("--" + curr_ins_dna)
            if curr_type == "y3":
                # 3n deletion
                curr_del_pro = self.dt.decodeAA(pro_num[curr_i-1])
                curr_ins_dna = ""
                align_pro.append("&"+curr_del_pro+"&")
                align_dna.append("---")
            if curr_type == "z1":
                # 3n - 2 insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                if nx_type == "w":
                    curr_ins_dna = dna_seq[curr_j-2]
                    curr_del_pro = "-"
                if nx_type == "z1":
                    curr_ins_dna = dna_seq[curr_j-4:curr_j-1]
                    curr_del_pro = "---"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)            
            if curr_type == "z2":
                # 3n - 1 insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                if nx_type == "w":
                    curr_ins_dna = dna_seq[curr_j-3:curr_j-1]
                    curr_del_pro = "--"
                if nx_type == "z2":
                    curr_ins_dna = dna_seq[curr_j-4:curr_j-1]
                    curr_del_pro = "---"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)
            if curr_type == "z3":
                # 3n insertion
                nx_i = curr_tracer.x
                nx_type = curr_tracer.prevMat
                nx_j = curr_tracer.y
                curr_ins_dna = dna_seq[curr_j-4:curr_j-1]
                curr_del_pro = "---"
                align_pro.append(curr_del_pro)
                align_dna.append(curr_ins_dna)
            nx_type = curr_tracer.prevMat
            
            if nx_type == None:
                break
            curr_type = nx_type
            curr_i = curr_tracer.x
            curr_j = curr_tracer.y
            #print(nx_type, curr_i, curr_j)
            if nx_type == "w":
                curr_tracer = w[curr_i][curr_j]
            if nx_type == "X":
                curr_tracer = X[curr_i][curr_j]
            if nx_type == "y1":
                curr_tracer = y1[curr_i][curr_j]
            if nx_type == "y2":
                curr_tracer = y2[curr_i][curr_j]
            if nx_type == "y3":
                curr_tracer = y3[curr_i][curr_j]
            if nx_type == "z1":
                curr_tracer = z1[curr_i][curr_j]
            if nx_type == "z2":
                curr_tracer = z2[curr_i][curr_j]
            if nx_type == "z3":
                curr_tracer = z3[curr_i][curr_j]
            if nx_type == "k":
                curr_tracer = k[curr_i][curr_j]
            if nx_type == "Y":
                curr_tracer = Y[curr_i][curr_j]
        head_pro = pro_seq[:curr_i-1]
        head_dna = dna_seq[:curr_j-1]
        return head_pro , head_dna , align_pro, align_dna, tail_pro, tail_dna, res


p = Parameters()
pg_parameters_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\data_analysis\final.txt"
with open(pg_parameters_file, "r") as f:
    all_lines = f.readlines()
    p.setTransitions(all_lines[-4])
    p.setPhi(all_lines[-3])
    p.setPsi(all_lines[-2])
    p.setPi(all_lines[-1])
    p.reNormalize()
    p.getScore()
    # p.fixTest()
rd_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\5000random_seq.txt"
rd.seed(10)
dt = DataTool()
with open(rd_file, "w+") as f:
    for i in range(5000):
        dna_seq_head = ">rd_dna_" + str(i) + "\n"
        pro_seq_head = ">rd_pro_" + str(i) + "\n"
        rd_pro = dt.genRandPro(200)
        rd_dna = dt.genRandDNA(500)
        f.write(pro_seq_head)
        f.write(rd_pro+"\n")
        f.write(dna_seq_head)
        f.write(rd_dna+"\n")

exit(0)
test_pro = "MSFPKYKPSSLRTLPETLDPAEYNISPETRRAQAERLAIRAQLKREYLLQYNDPNRRGLIENPALLRWAYARTINVYPNFRPTPKNSLMGALCGFGPLIFIYYIIKTERDRKEKLIQEGKLDRTFHLSYGTCGTTCCCAAAGTATACTC"
test_dna = "CTTCACCCCTGGCCACTATGCCCTCAAACCAGCCGAATAGGACAGGTCGCCAAGATTGTCTTCGGAAACCCAGCAGGCGCAAGCCAAGTGGTTGACCATTAAGAGCCTGGCTGAAACCTGCTTACCCTAACTGCCGAGGGCTCATCGAAAATCCTGCCTTGATTCCTTGGACCTATGCAAGATCAGCAAGTGTCTTTCCTAATTTCAGACCCACTCCTAAAAACTCACTCTTGGGAGCTCTGTGTGGAATTGGGCCCCTCTTCTTCTGGTATTATGTTTTCAAAACTGACAGGGATAGGAAAGAAAAACTTATCTGGGAAGGAAAATTGGGTCCAACATTTAACCTGTCATATTAAGTCTGG"

tp = "MYATMLWDQL"
td = "atgctatacggctctggatcagctta".upper()
#t = p.alignScore2(test_pro, test_dna)
# p.fixTest()
s = p.alignPattern(test_pro, test_dna)
print(s)

#test_dna = test_dna[::-1]
seq_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10_test.txt"

#exit(0)

dt = DataTool()
rd.seed(10)
with open(r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\gen_score1.txt", "w") as f:
    for i in range(5000):
        rd_pro = dt.genRandPro(200)
        rd_dna = dt.genRandDNA(500)
        k = p.alignScore(rd_pro, rd_dna)
        print(i, k)
        f.write(str(k)+"\n")

#with open(seq_file, "r") as f:
#    all_seqs = f.readlines()
#    tot = len(all_seqs)
#    with open(r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\test_score1.txt", "w") as g:
#        for i in range(tot//4):
#            pro_seq = all_seqs[4*i+1].rstrip()
#            dna_seq = all_seqs[4*i+3].rstrip().upper()
#            #s = p.alignPattern2(pro_seq, dna_seq)
#           t = p.alignScore(pro_seq, dna_seq)
#            u = p.alignScore(pro_seq[::-1], dna_seq)
#            # g.write(str(u) + " "+str(t)+"\n")
#            print(i, u, t)