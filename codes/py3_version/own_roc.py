import numpy as np
import matplotlib.pyplot as plt

def metrics(data, y, trueY):
    # 预测结果
    predicts = data[:] >= y
    print(predicts, trueY)
    # true pos
    tp = sum([(predicts[x] == 1 and trueY[x] == 1) for x in range(len(data))])
    # false pos
    fp = sum([( predicts[x] == 1 and trueY[x] == 0) for x in range(len(data))])
    # false neg
    fn = sum([( predicts[x] == 0 and trueY[x] == 1) for x in range(len(data))])
    # true neg
    tn = sum([( predicts[x] == 0 and trueY[x] == 0) for x in range(len(data))])
    print(tp, fp, fn, tn, y)
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
        params = xy(metrics(data, y, trueY))
        XY[i, 0] = params['fpr']
        XY[i, 1] = params['tpr']
        print(i, params)
    plt.plot(XY[:, 0], XY[:, 1], label=label)
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])

test_res = np.array([1, 1, 0, 0, 0])
real = np.array([1, 1, 0, 0, 0])
draw(test_res, real, [0, 0.2, 0.4, 0.6, 0.8, 1.0], "test")
plt.legend()
plt.show()
    