import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics 

import numpy as np
import matplotlib.pyplot as plt

def Metrics(data, y, trueY):
    # 预测结果
    predicts = data[:] >= y
    #print(predicts, trueY)
    # true pos
    tp = sum([(predicts[x] == 1 and trueY[x] == 1) for x in range(len(data))])
    # false pos
    fp = sum([( predicts[x] == 1 and trueY[x] == 0) for x in range(len(data))])
    # false neg
    fn = sum([( predicts[x] == 0 and trueY[x] == 1) for x in range(len(data))])
    # true neg
    tn = sum([( predicts[x] == 0 and trueY[x] == 0) for x in range(len(data))])
    #print(tp, fp, fn, tn, y)
    return {
        'tp':tp,
        'fp':fp,
        'tn':tn,
        'fn':fn
    }

def xy(params):
    fpr = params['fp'] / (params['fp']+params['tn'])
    tpr = params['tp'] / (params['tp']+params['fn'])
    return {
        'fpr':fpr,
        'tpr':tpr
    }


def draw(data, trueY, thresholds, label):
    N = len(thresholds)
    XY = np.zeros((N, 2))
    for i in range(N):
        y = thresholds[i]
        params = xy(Metrics(data, y, trueY))
        XY[i, 0] = params['fpr']
        XY[i, 1] = params['tpr']
        #print(i, params)
    plt.plot(XY[:, 0], XY[:, 1], label=label)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])

def get_score(num, rd_scores):
    tot = len(rd_scores)
    l, r = 0, tot - 1
    idx = -1
    while l <= r:
        mid = (l + r) // 2
        if rd_scores[mid] >= num:
            idx = mid
            r = mid - 1
        else:
            l = mid + 1
    if idx == -1:
        return 1.0
    else:
        return float(idx) / float(tot)

rd_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\gen_score1.txt"
test_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\codes\py3_version\test_score1.txt"
rd_fasty_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\rd_output.file"
pos_fasty_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\test_output.file"
neg_fasty_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\n_test_output.file"
rd_blastx_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\rd_output.file 2"
neg_blastx_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\rtest_output.file"
pos_blastx_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\test_output.file 2"
test_seq_file = r"C:\Users\InYuo\Documents\GitHub\protein-dna-align-EM\data\pg_500_hg19_presuf10_test.txt"
lengthSTAT = {}
tot_list = []
with open(test_seq_file, "r") as f:
    all_seqs = f.readlines()
    tot = len(all_seqs)
    for i in range(tot//4):
        pro_seq = all_seqs[4*i+1].rstrip()
        pro_hd = all_seqs[4*i].rstrip()
        dna_hd = all_seqs[4*i+2].rstrip()
        dna_seq = all_seqs[4*i+3].rstrip().upper()
        new_head = pro_hd+'_'+dna_hd
        lengthSTAT[new_head] = len(dna_seq)*len(pro_seq)
        tot_list.append(new_head)

new_stat = dict(sorted(lengthSTAT.items(), key=lambda item: item[1]))

idx = 0
short_part = {}
for key, length in new_stat.items():
    short_part[key] = length
    idx += 1
    if idx >= 50:
        break

def get_fasty_scores(fasty_file):
    curr_dict = {}
    with open(fasty_file, "r") as f:
        all_data = f.readlines()
        for line in all_data:
            curr_list = line.rstrip().split("\t")
            curr_score = float(curr_list[-1])
            curr_name = curr_list[0] + "_" + curr_list[1]
            if curr_name not in curr_dict:
                curr_dict[curr_name] = curr_score
            else:
                curr_dict[curr_name] = max(curr_dict[curr_name], curr_score)
    res = [y for x, y in curr_dict.items()]
    return res, curr_dict

def get_half_fasty_scores(fasty_file):
    with open(fasty_file, "r") as f:
        all_data = f.readlines()
        curr_dict = {}
        for line in all_data:
            curr_list = line.rstrip().split("\t")
            curr_score = float(curr_list[-1])
            curr_name = ">"+curr_list[1].rstrip() + "_>" + curr_list[0].rstrip()
            print("#", curr_name)
            if curr_name in short_part:
                if curr_name not in curr_dict:
                    curr_dict[curr_name] = curr_score
                else:
                    curr_dict[curr_name] = max(curr_dict[curr_name], curr_score)
    res = [y for x, y in curr_dict.items()]
    return res, curr_dict


def get_blastx_scores(blastx_file, keys):
    curr_dict = {}
    with open(blastx_file, "r") as f:
        all_data = f.readlines()
        for line in all_data:
            curr_list = line.rstrip().split("\t")
            curr_score = float(curr_list[-1])
            curr_name = curr_list[0] + "_" + curr_list[1]
            if curr_name not in curr_dict:
                curr_dict[curr_name] = curr_score
            else:
                curr_dict[curr_name] = max(curr_dict[curr_name], curr_score)
    for key in keys:
        if key not in curr_dict.keys():
            curr_dict[key] = 0
            print(key)
    res = [y for x, y in curr_dict.items()]
    return res

def get_half_blastx_scores(blastx_file, keys):
    curr_dict = {}
    with open(blastx_file, "r") as f:
        all_data = f.readlines()
        for line in all_data:
            curr_list = line.rstrip().split("\t")
            curr_score = float(curr_list[-1])
            curr_name = ">"+curr_list[1] + "_>" + curr_list[0]
            if curr_name in keys:
                if curr_name not in curr_dict:
                    curr_dict[curr_name] = curr_score
                else:
                    curr_dict[curr_name] = max(curr_dict[curr_name], curr_score)
    for key in keys:
        if key not in curr_dict.keys():
            curr_dict[key] = 0
            print(key)
    res = [y for x, y in curr_dict.items()]
    return res



rd_fasty_score, all_rd_keys = get_fasty_scores(rd_fasty_file)
pos_fasty_score, all_pos_keys = get_fasty_scores(pos_fasty_file)
neg_fasty_score, all_neg_keys = get_fasty_scores(neg_fasty_file)
rd_fasty_score.sort()

pos_blastx_score = get_blastx_scores(pos_blastx_file, all_pos_keys)
rd_blastx_score = get_blastx_scores(rd_blastx_file, all_rd_keys)
neg_blastx_score = get_blastx_scores(neg_blastx_file, all_neg_keys)
rd_blastx_score.sort()
all_blast_scores = [get_score(x, rd_blastx_score) for x in neg_blastx_score]
szz = len(pos_blastx_score)
szz2 = len(neg_blastx_score)
print(szz, szz2)

for i in range(szz):
    all_blast_scores.append(get_score(pos_blastx_score[i], rd_blastx_score))


sz = len(pos_fasty_score)
sz2 = len(neg_fasty_score)
all_scores = [get_score(x, rd_fasty_score) for x in neg_fasty_score]
for i in range(sz):
    all_scores.append(get_score(pos_fasty_score[i], rd_fasty_score))

rscores = []
with open(rd_file, "r") as f:
    rscores = [float(x.rstrip()) for x in f.readlines()]
rscores.sort()



test_scores = []
test_reverse = []
with open(test_file, "r") as f:
    for p in f.readlines():
        curr_t = float(p.rstrip().split()[2])
        curr_r = float(p.rstrip().split()[1])
        test_scores.append(curr_t)
        test_reverse.append(curr_r)

test_half_scores = []
test_half_reverse = []
with open(test_file, "r") as f:
    for p in f.readlines():
        index = int(p.rstrip().split()[0])
        curr_t = float(p.rstrip().split()[2])
        curr_r = float(p.rstrip().split()[1])
        curr_name = tot_list[index]
        if curr_name in short_part:
            test_half_scores.append(curr_t)
            test_half_reverse.append(curr_r)


fasty_half, fasty_dict = get_half_fasty_scores(pos_fasty_file)
blastx_half = get_half_blastx_scores(pos_blastx_file, fasty_dict)
print(blastx_half)
blastx_half_r = get_half_blastx_scores(neg_blastx_file, fasty_dict)
fasty_half_r, _ = get_half_fasty_scores(neg_fasty_file)

print(len(test_half_scores), len(test_half_reverse))
print(len(blastx_half), len(blastx_half_r))
print(len(fasty_half), len(fasty_half_r))
test_half = test_half_scores + test_half_reverse
blastx_half = blastx_half + blastx_half_r
fasty_half = fasty_half + fasty_half_r
true_half = [1 for x in range(50)] + [0 for x in range(50)]
print(blastx_half)
fpr1, tpr1, thresholds = metrics.roc_curve(true_half, np.array(test_half))
fpr2, tpr2, thresholds = metrics.roc_curve(true_half, np.array(fasty_half))
fpr3, tpr3, thresholds = metrics.roc_curve(true_half, np.array(blastx_half))
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr1, tpr1, label='our alignment algorithm')
plt.plot(fpr2, tpr2, label='fasty algorithm')
plt.plot(fpr3, tpr3, label = 'blastx algorithm')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
print(metrics.auc(fpr1, tpr1))
print(metrics.auc(fpr2, tpr2))
print(metrics.auc(fpr3, tpr3))
plt.legend(loc='best')
plt.savefig("50_raw_roc.png", bbox_inches='tight')
plt.show()




exit(0)

size = len(test_scores)
trueY = [0 for i in range(size)]

scores = [get_score(x, rscores) for x in test_reverse]
for i in range(size):
    trueY.append(1)
    t = get_score(test_scores[i], rscores)
    #print(i, t)
    scores.append(t)
#print(max(pos_blastx_score), max(pos_fasty_score), max(test_scores))
not_compressed_blastx = neg_blastx_score + pos_blastx_score
not_compressed_fasty = neg_fasty_score + pos_fasty_score
not_compressed = test_reverse + test_scores
#thresholds = np.array([x-0.99 for x in range(552)])
#l1 = draw(not_compressed, trueY, thresholds, "our algorithm")
#l2 = draw(not_compressed_blastx, trueY, thresholds, "blastx")
#l3 = draw(not_compressed_fasty, trueY, thresholds, "fasty")
#thresholds = np.array([(1-0.002*(x-1)) for x in range(1,2)])
#l1 = draw(scores, trueY, thresholds, "our algorithm")
#l2 = draw(all_blast_scores, trueY, thresholds, "blastx")
#l3 = draw(all_scores, trueY, thresholds, "fasty")

#plt.legend(loc="best")
#plt.show()
#exit(0)

print(pos_blastx_score)
print(test_scores)
print(len(test_scores), len(test_reverse))
y = np.array(trueY)
tscores = np.array(scores)

fpr2, tpr2, thresholds = metrics.roc_curve(y, np.array(not_compressed_fasty))

fpr3, tpr3, thresholds = metrics.roc_curve(y, np.array(not_compressed_blastx))
fpr, tpr, thresholds = metrics.roc_curve(y, np.array(not_compressed))
#thresholds = np.array([x for x in range(551)])
#print(thresholds)
#plt.figure(1)
#metrics.roc_curve()
#ffpr, ttpr, tthresholds = me.roc_curve(y, scores)
#print(tthresholds)
#exit(0)
#ffpr2, ttpr2, tthresholds = metrics.roc_curve(y, np.array(not_compressed_fasty))
#ffpr3, ttpr3, tthresholds = metrics.roc_curve(y, np.array(not_compressed_blastx))
#tthresholds = np.array([(1-0.0002*(x-1)) for x in range(5001)])
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr, label='our alignment algorithm')
plt.plot(fpr2, tpr2, label='fasty algorithm')
plt.plot(fpr3, tpr3, label = 'blastx algorithm')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
print(metrics.auc(fpr, tpr))
print(metrics.auc(fpr2, tpr2))
print(metrics.auc(fpr3, tpr3))
plt.legend(loc='best')
plt.savefig("full_raw_roc.pdf", bbox_inches='tight')

plt.show()
