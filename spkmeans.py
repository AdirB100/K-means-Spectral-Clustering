import numpy as np
import sys
import pandas as pd
import mykmeanssp as km

np.random.seed(0)

## vectorsList should be np.vector
def kmeanspp(vectorsList, k):
    i = 0
    N = len(vectorsList)
    D = np.array([0.0 for q in range(N)])
    P = np.array([0.0 for q in range(N)])
    indexList = []
    centroids = np.array([[0.0 for n in range(len(vectorsList[0]))] for m in range(k)])
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
        if i == k - 1:
            break
    return centroids, indexList

def ZeroCheck(x):
    if x<0 and x>-0.00005:
        return 0
    return x

### condition value for epsilon
try:
    if len(sys.argv)!=4:
        raise ValueError
    k = float(sys.argv[1])
    if k <0 or int(k)!=k or k==1:
        raise ValueError
    goal_options=["spk","wam","ddg","lnorm","jacobi"]
    goal=sys.argv[2]
    if goal not in goal_options:
        raise ValueError
    file_name=sys.argv[3]
    df = pd.read_csv(file_name, header=None)

except:
    print("Invalid Input!")
    sys.exit(1)

vectorsList = df.to_numpy()
d = len(vectorsList[0])
N = len(vectorsList)
k = int(k)
if k >= N: ## checks if k is valid
    print("Invalid Input!")
    exit(1)

if goal == "spk":
    T_mat=np.array(km.getTmatPy(N, d, k, vectorsList.tolist()))
    k=len(T_mat[0])
    if k<=1:
        print("An Error Has Occurred")
        sys.exit(1)
    # finding first centroids and indices by kmeans++
    centroids, indices = kmeanspp(T_mat, k)

    ### This loop prints the indices of the centroids
    for i in range(len(indices)):
        if i == len(indices) - 1:
            print(indices[i])
            break
        print(str(indices[i]) + ',', end='')

    final_centroids=km.fit(centroids.tolist())
    if final_centroids==None:
        print("An Error Has Occurred")
        sys.exit(1)
    for i in range(len(final_centroids)):
        for j in range(len(final_centroids[0])):
            if j == len(final_centroids[0]) - 1:
                print("{:.4f}".format(final_centroids[i][j]))
                break
            print("{:.4f}".format(final_centroids[i][j]) + ',', end='')
   
else:
    mat=km.notSpk(N, d, k, goal, vectorsList.tolist())
    if mat==None:
        print("An Error Has Occurred")
        sys.exit(1)
    if goal=="jacobi":
        eigenVals=mat.pop(0)
        for i in range(len(eigenVals)):
            if i != len(eigenVals) - 1:
                print("{:.4f}".format(ZeroCheck(eigenVals[i]))+',',end='')
            else:
                print("{:.4f}".format(ZeroCheck(eigenVals[i])))
    for i in range(len(mat)):
        for j in range(len(mat[0])):
            if j == len(mat[0]) - 1:
                print("{:.4f}".format(mat[i][j]))
                break
            print("{:.4f}".format(mat[i][j]) + ',', end='')

sys.exit(0)

