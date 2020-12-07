#coding=utf-8

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import log, exp

# pg_parameters_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\parameter_log_SMALL_PG100.txt"
pg_parameters_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\parameter_log.txt"
# pg_parameters_file = r'/Users/yinyao/mt/protein-dna-align-EM/codes/cpp_version/parameter_log_pg_nopt.txt'
# pg_parameters_file = r'/Users/yinyao/mt/protein-dna-align-EM/codes/cpp_version/parameter_log.txt'
# pg_error_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\error_log_SMALL_PG100.txt"
# pg_error_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\error_log.txt"
pg_data_df = pd.DataFrame(columns=["epoch", "prob", "updateMethod", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"])
#pg_error_df = pd.DataFrame(columns=["epoch", "prob", "type", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"])

pg_line_map = {}
pg_error_map = {}
pg_line_key = ["epoch", "prob", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "gamma", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"]

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
    vals = np.round(np.array(vals), 5)
    R, C = len(row), len(col)
    idx = pd.Index(row)
    df = pd.DataFrame(np.random.randn(R, C), index=idx, columns=col)
    
    # 根据行数列数设置表格大小
    figC, figR = 2.25*C, R
    fig = plt.figure(figsize=(figC, figR))
    
    # 设置fig并去掉边框
    ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    the_table=plt.table(cellText=vals, rowLabels=df.index, colLabels=df.columns, colWidths=[0.005]*vals.shape[1], rowLoc='center', loc='center',cellLoc='center')
    the_table.set_fontsize(15)
    
    # 伸缩表格大小常数
    the_table.scale(figR/R*2, figC/C*1.5)
    plt.savefig("testdoc.jpg")



'''
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
'''

"""
This is the 0 epoch: 
cnts: 1079789.00714693	21387.4396065484	98623.7001027322	6629	1091441.8814422	20571.9850908884	14700.4921535838	18.9191545538601	152.403784852831	13137.1081463581	1410.98022237288	14240.7284560161	171.322939406691	99.0333659022995	1049.77450037673	11190.4519855665	12339.2598518455	1159.20033136252	109.425830985798	10.3924650834985	
Overall: -148014.366331961
0.759894445124435	0.993930836715862	0.745497183636171	0.0932727619254954	0.111121330277092	0.745497183636171	0.092787979548694	0.0829322891004948	0.10101584601602	0.105630808440827	0.090156076639369	0.118496949446549
0.0442408874577312	0.0763118312089171	0.0574485466210105	0.0421403032197225	0.00353279328008127	0.0582099790280728	0.117376462619846	0.031323318574245	0.0128623456777058	0.081843025397497	0.0735904788480046	0.122404023135783	0.0307741447468958	0.0304267685273387	0.0395404531659081	0.0629938275795671	0.0607642523861573	0.00566940843785488	0.0165091769823467	0.029013186835476	0.00302478626983846	
0.271147458326026	0.14847203321557	0.405973698665912	0.174406809792491	
blabla
"""


def strToCnts(s):
    cnts_mp = {}
    cnts_keys = ["J_d.cnt", "J_i.cnt", "M.cnt", "A.cnt", "K_d.cnt", "K_i.cnt", 
        "F_d.cnt", "X_d.cnt", "B_d.cnt", "D_d.cnt", "E_d.cnt", "G_d.cnt", "H_d.cnt",
        "B_i.cnt", "E_i.cnt", "D_i.cnt", "F_i.cnt", "G_i.cnt", "H_i.cnt", "X_i.cnt"]
    curr_list = s.rstrip().split("\t")
    # print(curr_list)
    for i in range(len(curr_list)):
        cnts_mp[cnts_keys[i]] = float(curr_list[i])
    return cnts_mp


def strToPara(s):
    para_list = s.rstrip().split("\t")
    para_key_list = ["omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "gamma", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d"]
    para_mp = {}
    for i in range(len(para_key_list)):
        para_mp[para_key_list[i]] = float(para_list[i])
    return para_mp


def strToPi(s):
    pass


def strToPsi(s):
    pass


def strToPhi(s):
    pass


def reCalulate(cnts, omega_i, omega_d, prob):
    # (B_i.cnt + D_i.cnt + E_i.cnt)*log_alpha_i
    cnts_mp = strToCnts(cnts)
    prob += log(omega_i)*(cnts_mp["J_i.cnt"]+cnts_mp["K_i.cnt"])+log(omega_d)*(cnts_mp["J_d.cnt"] + cnts_mp["K_d.cnt"])
    prob += 2*(cnts_mp["A.cnt"]*log(1-omega_i) + cnts_mp["A.cnt"]*log(1-omega_d)) 
    return prob


def handle(mmap):
    newmp = {}
    for key in ["epoch", "prob", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"]:
        newmp[key] = [mmap[key]]
    df2 = pd.DataFrame.from_dict(newmp)
    return df2


# load data file
with open(pg_parameters_file, "r") as f:
    # read 7 lines
    all_lines = f.readlines()
    curr_list = []
    for i in range(len(all_lines)):
        curr = all_lines[i].rstrip()
        if i % 10 == 0:
            # epoch number
            ep = int(curr.split(" ")[3])
            pg_line_map["epoch"] = ep
        elif i % 10 == 1:
            cnts = curr[5:]
            pg_line_map["cnts"] = cnts
        elif i % 10 == 5:
            prob = float(curr[9:])
            pg_line_map["prob"] = prob
        elif i % 10 == 6:
            curr_para_map = strToPara(curr)
            for k, v in curr_para_map.items():
                pg_line_map[k] = v
        elif i % 10 == 7:
            pg_line_map["phi"] = curr
        elif i % 10 == 8:
            pg_line_map["psi"] = curr
        elif i % 10 == 9:
            pg_line_map["pi"] = curr
            #pg_line_map["prob"] = reCalulate(pg_line_map["cnts"], pg_line_map["omega_i"],pg_line_map['omega_d'], pg_line_map["prob"])
            pg_data_df = pg_data_df.append(pg_line_map, ignore_index=True)

"""
can not find optimal deletion parameters: 
cnts: 9546.21619881218	11416.4153027483	127406.692275202	1182.01	9371.78031208009	13377.6857610754	13572.3512126715	0.0100012342487303	1107.35920518579	12459.5445238472	5.46748363857974	13084.9781391371	1107.35920642003	507.929986453203	696.082188257924	14286.0884387362	15490.0806134473	1205.78977625312	509.717587995193	1.79760154198995	
Overall: -139048.236760013
0.912953477361546	0.888919360682289	0.808155682959386	0.0982555360069311	0.0860911521277267	0.808155682959386	0.0748182067392156	0.00244717908272717	0.458536873206279	0.0783789504655535	1.00151367731258e-009	0.994973771113063
0.099800713754121	0.0546631313727216	0.0270001546551439	0.0394640430965702	0.0138347487047729	0.0398933038065695	0.0683476404580205	0.0504439127935746	0.0137587985920525	0.03013490519776	0.0708836203966992	0.0767093372206357	0.083812964991124	0.0323511237039265	0.0585818861989649	0.0726473210975572	0.0476597793325595	0.00951974703937268	0.0267082364379659	0.0474700423262555	0.0363145888236322	
0.323491356296386	0.175012882209597	0.31136515553601	0.190130605958007	
blabla
"""
"""
with open(pg_error_file, "r") as f:
    # read 7 lines
    all_lines = f.readlines()
    curr_list = []
    for i in range(len(all_lines)):
        curr = all_lines[i].rstrip()
        if i % 7 == 0:
            # deletion/insertions
            ttype = curr.split(" ")[4]
            pg_error_map["type"] = ttype
        elif i % 7 == 3:
            curr_error_map = strToPara(curr)
            for k, v in curr_error_map.items():
                pg_error_map[k] = v
        elif i % 7 == 4:
            pg_error_map["phi"] = curr
        elif i % 7 == 5:
            pg_error_map["psi"] = curr
        elif i % 7 == 6:
            pg_error_df = pg_error_df.append(pg_error_map, ignore_index = True)

# have pg_error_df, pg_data_df, need to have some label

def check(gamma, error_df):
    curr_res = ""
    for idx, row in error_df.iterrows():
        if row["gamma"] == gamma:
            curr_res += row["type"]
    return curr_res

for idx, row in pg_data_df.iterrows():
    row["updateMethod"] = check(row["gamma"], pg_error_df)
    #print(idx, row["updateMethod"])


"""


def get_fig(df, epoch_start, epoch_end, col, name):
    x1 = [i for i in range(epoch_start, epoch_end + 1, 1)]
    y1 = []
    for idx, row in df.loc[epoch_start:epoch_end].iterrows():
        y1.append(log(row[col]))
    l1 = plt.plot(x1, y1, 'r--', label ='type1')
    plt.legend()
    plt.show()

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

p = Parameters()

with open(pg_parameters_file, "r") as f:
    all_lines = f.readlines()
    p.setTransitions(all_lines[-4])
    p.setPhi(all_lines[-3])
    p.setPsi(all_lines[-2])
    p.setPi(all_lines[-1])
    p.reNormalize()
    p.getScore()
    # print(p.s[0][1])

#for idx, row in pg_data_df.iterrows():
    #print(idx, row['prob'])

get_fig(pg_data_df, 0, 349, "alpha_d", "")
dt = DataTool()
row = dt.aaList
col = dt.tripletList

for i in range(64):
    curr_t = dt.decodeTriplet(i)
    curr_v = -10.0
    curr_idx = -1
    for j in range(21):
        if p.s[j][i] > curr_v:
            curr_v = p.s[j][i]
            curr_idx = j
    curr_aa = dt.decodeAA(curr_idx)
    if curr_t not in dt.codonTable[curr_aa]:
        print("Not match ", curr_t, curr_aa)
    else:
        print("Match ", curr_t, curr_aa)

print("hello\n")
for i in range(21):
    curr_aa = dt.decodeAA(i)
    curr_v = -10.0
    curr_idx = -1
    for j in range(64):
        if p.s[i][j] > curr_v:
            curr_v = p.s[i][j]
            curr_idx = j
    curr_t = dt.decodeTriplet(curr_idx)
    if curr_t not in dt.codonTable[curr_aa]:
        print("Not match ", curr_t, curr_aa)
    else:
        print("Match ", curr_t, curr_aa)
    # print(curr_t, dt.decodeAA(curr_idx))
# plot_table(row, col, p.s)

def patternCollect(strList):
    stat = [0 for i in range(64)]
    dt = DataTool()
    for s in strList:
        for i in range(len(s)-2):
            stat[dt.encodeTriplet(s[i:i+3])] += 1
    return [ x/sum(stat) for x in stat]

seqFile = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\small_test_pg.txt"

with open(seqFile, "r") as f:
    tot = f.readlines()
    curr = []
    for i in range(len(tot)):
        if i % 4 == 3 and not tot[i].startswith("N"):
            curr.append(tot[i].rstrip().upper())
    stat = patternCollect(curr)
    p.getNewScore(stat)
    print("###")
    for i in range(64):
        curr_t = dt.decodeTriplet(i)
        curr_v = -10.0
        curr_idx = -1
        for j in range(21):
            if p.s[j][i] > curr_v:
                curr_v = p.s[j][i]
                curr_idx = j
        curr_aa = dt.decodeAA(curr_idx)
        if curr_t not in dt.codonTable[curr_aa]:
            print("Not match ", curr_t, curr_aa)
        else:
            print("Match ", curr_t, curr_aa)
    for i in range(21):
        curr_aa = dt.decodeAA(i)
        curr_v = -10.0
        curr_idx = -1
        for j in range(64):
            if p.s[i][j] > curr_v:
                curr_v = p.s[i][j]
                curr_idx = j
        curr_t = dt.decodeTriplet(curr_idx)
        if curr_t not in dt.codonTable[curr_aa]:
            print("Not match ", curr_t, curr_aa)
        else:
            print("Match ", curr_t, curr_aa)
    plot_table(row, col, p.s)
    