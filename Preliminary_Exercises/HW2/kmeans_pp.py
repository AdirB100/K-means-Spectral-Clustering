import numpy as np
import sys
import pandas as pd
import mykmeanssp as mk

np.random.seed(0)

## vectorsList should be np.vector
def kmeanspp(vectorsList, K):
    i = 0
    N = len(vectorsList)
    D = np.array([0.0 for q in range(N)])
    P = np.array([0.0 for q in range(N)])
    indexList = []
    centroids = np.array([[0.0 for n in range(len(vectorsList[0]))] for m in range(K)])
    randIndex = np.random.choice(np.arange(len(vectorsList)))
    indexList.append(randIndex)
    centroids[0] = vectorsList[randIndex]
    while (True):
        for l in range(N):
            D[l] = np.inf
            for j in range(i + 1):
                D[l] = np.minimum(D[l], (np.linalg.norm(vectorsList[l] - centroids[j])) ** 2)
        for l in range(N):
            P[l] = D[l] / (np.sum(D))
        i += 1
        randIndex = np.random.choice(np.arange(len(vectorsList)), p=P)
        indexList.append(randIndex)
        centroids[i] = vectorsList[randIndex]
        if i == K - 1:
            break
    return centroids, indexList


### condition value for epsilon
try:
    K = int(sys.argv[1])
    if K <= 1:
        raise ValueError
    if len(sys.argv) < 6:
        max_iter = 300
        eps = float(sys.argv[2])
        if eps < 0:
            raise ValueError
        file_name_1 = sys.argv[3]
        file_name_2 = sys.argv[4]
    else:
        max_iter = int(sys.argv[2])
        eps = float(sys.argv[3])
        if eps < 0:
            raise ValueError
        if max_iter <= 0:
            raise ValueError
        file_name_1 = sys.argv[4]
        file_name_2 = sys.argv[5]
    df1 = pd.read_csv(file_name_1, header=None)
    df2 = pd.read_csv(file_name_2, header=None)
except:
    print("Invalid Input!")
    sys.exit(1)
df = pd.merge(df1, df2, how='inner', on=0, sort=True)
df.drop(0, axis=1, inplace=True)
vectorsList = np.array(df.values.tolist())
centroids, indices = kmeanspp(vectorsList, K)
d = len(vectorsList[0])
N = len(vectorsList)
centroids = mk.fit(K, N, d, max_iter, eps, vectorsList.tolist(), centroids.tolist())
if centroids == None:
    print('An Error Has Occurred')
for i in range(len(indices)):
    if i == len(indices) - 1:
        print(indices[i])
        break
    print(str(indices[i]) + ',', end='')
for i in range(K):
    for j in range(d):
        if j == d - 1:
            print("{:.4f}".format(centroids[i][j]))
            break
        print("{:.4f}".format(centroids[i][j]) + ',', end='')
