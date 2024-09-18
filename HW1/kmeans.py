import sys


def kmeans(K, inputfile, maxiter):
    try:
        f = open(inputfile, 'r')
    except FileNotFoundError:
        print("Invalid Input!")
        sys.exit(1)
    except:
        print("An Error Has Occurred")
        sys.exit(1)
    vectorslist = make_veclist(f)
    if K >= len(vectorslist):
        print("Invalid Input!")
        sys.exit(1)
    f.close()
    centroids = [vectorslist[i] for i in range(K)]
    cluster_mem = [0 for i in range(len(vectorslist))]
    itercnt = 0
    norm = 0.001
    while itercnt < maxiter:
        itercnt += 1
        for i in range(len(vectorslist)):
            cluster_mem[i] = getclosestmeanindex(vectorslist[i], centroids)
        newcentroids = centroids.copy()
        for k in range(K):
            newcentroids[k] = calcnewmean(k, cluster_mem, vectorslist)
        norm = calcnorm(newcentroids, centroids)
        centroids = newcentroids
        if norm:
            break
    return centroids


def calcnorm(newcent, oldcent):
    val = 0
    d = len(newcent[0])
    for i in range(len(newcent)):
        subvec = [newcent[i][j] - oldcent[i][j] for j in range(d)]
        for ele in subvec:
            val += ele ** 2
        if val ** 0.5 >= 0.001:
            return False
    return True


def calcnewmean(k, cluster_mem, vectorlist):
    '''the function calculates the means by the conected vectors '''
    d = len(vectorlist[0])
    N = len(vectorlist)
    groupsize = 0
    sumvec = [0 for i in range(d)]
    for i in range(N):
        if cluster_mem[i] == k:
            sumvec = [sumvec[j] + vectorlist[i][j] for j in range(d)]
            groupsize += 1
    return [sumvec[i] / groupsize for i in range(d)]


def getclosestmeanindex(vector, centroids):
    '''for each vector compute the index of the closest cluster'''
    d = len(vector)
    minval = 0
    minindex = 0
    sub0 = [vector[j] - centroids[0][j] for j in range(d)]
    for ele in sub0:
        minval += ele ** 2
    for i in range(1, len(centroids)):
        val = 0
        sub = [vector[j] - centroids[i][j] for j in range(d)]
        for ele in sub:
            val += ele ** 2
        if (val < minval):
            minindex = i
            minval = val
    return minindex


def make_veclist(f):
    '''make vector list of floats from the file $ret=array of arrays, each internal array is a vector'''
    inparr = f.read().split('\n')
    for i in range(len(inparr) - 1, -1, -1):
        if len(inparr[i]) == 0:
            inparr.pop(i)
    vectorslist = [inparr[i].split(',') for i in range(len(inparr))]
    for i in range(len(vectorslist)):
        for j in range(len(vectorslist[i])):
            vectorslist[i][j] = float(vectorslist[i][j])
    return vectorslist


##### Main Code #####
try:
    K = int(sys.argv[1])
    if K <= 1:
        raise ValueError
    if len(sys.argv) < 5:
        maxiter = 200
        inputfile = sys.argv[2]
        outputfile = sys.argv[3]
    else:
        maxiter = int(sys.argv[2])
        if maxiter <= 0:
            raise ValueError
        inputfile = sys.argv[3]
        outputfile = sys.argv[4]
except:
    print("Invalid Input!")
    sys.exit(1)
centroids = kmeans(K, inputfile, maxiter)
try:
    f = open(outputfile, 'w')
except:
    print("An Error Has Occurred")
    sys.exit(1)
for centroid in centroids:
    for i in range(len(centroid) - 1):
        f.write("{:.4f}".format(centroid[i]) + ',')
    f.write("{:.4f}".format(centroid[len(centroid) - 1]))
    f.write('\n')
f.close()
sys.exit(0)
