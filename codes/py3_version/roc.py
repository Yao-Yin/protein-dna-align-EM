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
thresholds = np.array([(1-0.0002*(x-1)) if x != 0 else 2 for x in range(5001)])
print(thresholds)
plt.figure(1)
plt.plot([0, 1], [0, 1], 'k--')
plt.plot(fpr, tpr, label='Alignment algorithm')
plt.xlabel('False positive rate')
plt.ylabel('True positive rate')
plt.title('ROC curve')
print(metrics.auc(fpr, tpr))
plt.savefig("roc.pdf", bbox_inches='tight')
# plt.legend(loc='best')
# plt.show()
