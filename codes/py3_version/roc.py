import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import metrics

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
    return res

rd_fasty_score = get_fasty_scores(rd_fasty_file)
pos_fasty_score = get_fasty_scores(pos_fasty_file)
neg_fasty_score = get_fasty_scores(neg_fasty_file)
rd_fasty_score.sort()


print(get_fasty_scores(pos_fasty_file))
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

size = len(test_scores)
trueY = [0 for i in range(size)]

scores = [get_score(x, rscores) for x in test_reverse]
for i in range(size):
    trueY.append(1)
    t = get_score(test_scores[i], rscores)
    print(i, t)
    scores.append(t)

y = np.array(trueY)
tscores = np.array(scores)
fpr, tpr, thresholds = metrics.roc_curve(y, scores)
fpr2, tpr2, thresholds = metrics.roc_curve(y, np.array(all_scores))
thresholds = np.array([(1-0.0002*(x-1)) if x != 0 else 2 for x in range(5001)])
print(thresholds)
plt.figure(1)
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr, label='Alignment algorithm')
plt.plot(fpr2, tpr2, label='fasty algorithm')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
print(metrics.auc(fpr, tpr))
print(metrics.auc(fpr2, tpr2))
# plt.savefig("roc.pdf", bbox_inches='tight')
plt.legend(loc='best')
plt.show()
