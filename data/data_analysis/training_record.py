import matplotlib.pyplot as plt 
import pandas as pd 
from math import log
#pg_parameters_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\parameter_log_SMALL_PG100.txt"
pg_parameters_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\parameter_log.txt"
pg_error_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\error_log_SMALL_PG100.txt"
#pg_error_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\cpp_version\error_log.txt"
pg_data_df = pd.DataFrame(columns=["epoch", "prob", "updateMethod", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"])
pg_error_df = pd.DataFrame(columns=["epoch", "prob", "type", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"])

pg_line_map = {}
pg_error_map = {}
pg_line_key = ["epoch", "prob", "omega_i", "omega_d", "gamma", "alpha_i", "alpha_d", "gamma", "delta_i", "beta_i", "epsilon_i", "delta_d", "beta_d", "epsilon_d", "cnts", "psi", "phi", "pi"]

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
    #print(curr_list)
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

def reCalulate(cnts, alpha_i, prob):
    #(B_i.cnt + D_i.cnt + E_i.cnt)*log_alpha_i
    cnts_mp = strToCnts(cnts)
    prob += log(alpha_i)*(cnts_mp["B_i.cnt"]+cnts_mp["D_i.cnt"]+cnts_mp["E_i.cnt"])
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
        if i % 7 == 0:
            #epoch number
            ep = int(curr.split(" ")[3])
            pg_line_map["epoch"] = ep
        elif i % 7 == 1:
            cnts = curr[5:]
            pg_line_map["cnts"] = cnts
        elif i % 7 == 2:
            prob = float(curr[9:])
            pg_line_map["prob"] = prob
        elif i % 7 == 3:
            curr_para_map = strToPara(curr)
            for k, v in curr_para_map.items():
                pg_line_map[k] = v
        elif i % 7 == 4:
            pg_line_map["phi"] = curr
        elif i % 7 == 5:
            pg_line_map["psi"] = curr
        elif i % 7 == 6:
            pg_line_map["pi"] = curr
            #pg_line_map["prob"] = reCalulate(pg_line_map["cnts"], pg_line_map["alpha_i"], pg_line_map["prob"])
            pg_data_df = pg_data_df.append(pg_line_map, ignore_index = True)

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
        y1.append((row[col]))
    l1=plt.plot(x1,y1,'r--',label='type1')
    plt.legend()
    plt.show()

get_fig(pg_data_df, 0, 59, "beta_", "")