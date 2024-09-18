#include "spkmeans.h"

/*Calculate the Weighted Adjacencey matrix
Input: vectorsList, num of vectors, dimenstion of vectors
*/
double **calcWAM(double **vectorsList, int N, int d)
{
    double **wam;
    int i, j;
    wam = makeMatrixOfDoubles(N, N);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < i; j++)
        {
            wam[i][j] = calcWeight(vectorsList[i], vectorsList[j], d);
            wam[j][i] = wam[i][j];
        }
        wam[i][i] = 0.0;
    }
    return wam;
}

/* Calculate wij*/
double calcWeight(double *vec_a, double *vec_b, int d)
{
    int i;
    double val, *sub;
    sub = calloc(d, sizeof(double));
    assert(sub != NULL);
    val = 0.0;
    for (i = 0; i < d; i++)
    {
        sub[i] = vec_a[i] - vec_b[i];
        val += sub[i] * sub[i];
    }
    val = exp(-sqrt(val) / 2);
    free(sub);
    return val;
}

/* Calculate the Diagonal Degree Matrix power -1/2 (D^-0.5)
Input: Weighted Matrix and dimension of the matrix
*/
double **calcDDM(double **wam, int N)
{
    int i, z;
    double **ddm;
    ddm = makeMatrixOfDoubles(N, N);
    for (i = 0; i < N; i++)
    {
        for (z = 0; z < N; z++)
        {
            ddm[i][i] += wam[i][z];
        }
    }
    return ddm;
}

/*Calculate L_norm matrix*/
double **calcLnorm(double **ddm, double **wam, int N)
{
    int i, j;
    double **Lnorm;
    Lnorm = makeMatrixOfDoubles(N, N);
    diagMatMUL(ddm, wam, &Lnorm, N, 0);
    diagMatMUL(Lnorm, ddm, &Lnorm, N, 1);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            if (i == j)
                Lnorm[i][j] = 1 - Lnorm[i][j];
            else
                Lnorm[i][j] = -1 * Lnorm[i][j];
        }
    }
    return Lnorm;
}

/*matrix multipication inplace, where one of the matrix is diagonal.
 whereDiag=0 if mat1 (the left one) is diagonal else 1
 */
void diagMatMUL(double **mat1, double **mat2, double ***resMatPtr, int N, int whereDiag)
{
    int i, j;
    if (!whereDiag)
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; ++j)
            {
                (*resMatPtr)[i][j] = pow(mat1[i][i], -0.5) * mat2[i][j];
            }
        }
    }
    else
    {
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; ++j)
            {
                (*resMatPtr)[i][j] = mat1[i][j] * pow(mat2[j][j], -0.5);
            }
        }
    }
}

/* compute the eigenvalus and eigenvectors. A is Lnorm.
A matrix changes here */
double **jacobi(double **A, int N)
{
    int i, j, rotations;
    double **V, deltaOffDiag;
    PMat *currP;
    currP = (PMat *)malloc(sizeof(PMat));
    assert(currP != NULL);
    rotations = 0;
    V = makeMatrixOfDoubles(N, N);
    deltaOffDiag = 0;
    /* calculating the first value of offA^2- offA'^2 */
    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            deltaOffDiag += 2 * A[i][j] * A[i][j];
        }
    }
    for (i = 0; i < N; i++)
        V[i][i] = 1; /* V starts as a unit matrix*/
    while (rotations < MAXROTATIONS && deltaOffDiag > EPS)
    { /*This Loop diagonalize matrix A and simultaniously calculate V*/
        buildP(A, currP, N);
        deltaOffDiag = 2 * A[currP->i][currP->j] * A[currP->i][currP->j]; /*this is the only diffrence between the squared off diagonal elements*/
        TransformA(&A, currP, N);
        calcNewV(&V, currP, N);
        rotations++;
    }
    free(currP);
    return V;
}

void buildP(double **A, PMat *currP, int N)
{
    /* Check again later for edge cases. Especially when A is 1x1 or 0x0*/
    int i, j, l, m, sign;
    double t, theta, currMax;
    currMax = -1;
    for (l = 0; l < N; l++)
    { /* Finding i and j of max value in matrix A off diagonal*/
        for (m = l + 1; m < N; m++)
        {
            sign = A[l][m] >= 0 ? 1 : -1;
            if (currMax < (A[l][m] * sign))
            {
                currMax = (A[l][m] * sign);
                currP->i = l;
                currP->j = m;
            }
        }
    }
    i = currP->i;
    j = currP->j; /*for convinient*/
    theta = (A[j][j] - A[i][i]) / (2 * A[i][j]);
    sign = theta >= 0 ? 1 : -1;
    t = sign / (theta * sign + sqrt(theta * theta + 1));
    currP->c = 1 / sqrt(t * t + 1);
    currP->s = currP->c * t;
}

void TransformA(double ***APtr, PMat *currP, int N)
{
    int i, j, r;
    double c, s;
    double a_ri, a_rj, a_ii, a_jj, a_ij;
    i = currP->i;
    j = currP->j;
    c = currP->c;
    s = currP->s;
    a_ii = (*APtr)[i][i];
    a_jj = (*APtr)[j][j];
    a_ij = (*APtr)[i][j];
    for (r = 0; r < N; r++)
    {
        if (r != i && r != j)
        {
            a_ri = (*APtr)[r][i];
            a_rj = (*APtr)[r][j];
            (*APtr)[r][i] = c * a_ri - s * a_rj;
            (*APtr)[r][j] = c * a_rj + s * a_ri;
            (*APtr)[i][r] = (*APtr)[r][i];
            (*APtr)[j][r] = (*APtr)[r][j];
        }
    }
    (*APtr)[i][i] = c * c * a_ii + s * s * a_jj - 2 * s * c * a_ij;
    (*APtr)[j][j] = s * s * a_ii + c * c * a_jj + 2 * s * c * a_ij;
    (*APtr)[i][j] = (*APtr)[j][i] = 0;
}

void calcNewV(double ***VPtr, PMat *currP, int N)
{
    int i, j, m;
    double a_mi, a_mj, c, s;
    i = (currP->i);
    j = (currP->j);
    c = (currP->c);
    s = (currP->s);
    for (m = 0; m < N; m++)
    {
        a_mi = (*VPtr)[m][i];
        a_mj = (*VPtr)[m][j];
        (*VPtr)[m][i] = (c * a_mi) + ((-s) * a_mj);
        (*VPtr)[m][j] = (c * a_mj) + (s * a_mi);
    }
}

double **makeMatrixOfDoubles(int N, int d)
{
    int i;
    double *resSUB, **res;
    res = (double **)calloc(N, sizeof(double *));
    assert(res != NULL);
    resSUB = (double *)calloc(N * d, sizeof(double));
    assert(resSUB != NULL);
    for (i = 0; i < N; i++)
        res[i] = resSUB + i * d;
    return res;
}

int eigenGap(Eigen *eigenArr, int N)
{
    int i, k;
    double delta, maxDelta;
    k = 0;
    maxDelta = -1;
    for (i = 1; i <= (int)(N/2); i++)
    {
        delta = eigenArr[i].eigenVal - eigenArr[i-1].eigenVal;
        if (delta > maxDelta)
        {
            maxDelta = delta;
            k = i;
        }
    }
    return k;
}
/* Fill and sort Eigen Array*/
Eigen *eigenFill(double **Lnorm, double **V, int N)
{
    int i;
    double **VTransposed;
    Eigen *eigenArr;
    VTransposed = transpose(V, N);
    eigenArr = (Eigen *)malloc(N*sizeof(Eigen));
    assert(eigenArr != NULL);
    for (i = 0; i < N; i++)
    {
        eigenArr[i].eigenVal = Lnorm[i][i];
        eigenArr[i].eigenVec = VTransposed[i];
        eigenArr[i].orIndex = i; /*Can be deleted? */
    }
    qsort(eigenArr, N, sizeof(Eigen), comparator);
    return eigenArr;
}

/**used to sort EigenArray*/
int comparator(const void *p1, const void *p2)
{
    double delta;
    Eigen *eigen1, *eigen2;
    eigen1 = (Eigen *)p1, eigen2 = (Eigen *)p2;
    delta = ((eigen1->eigenVal) - (eigen2->eigenVal));
    if (delta == 0)
        return 0;
    return (delta > 0 ? 1 : -1);
}

double **transpose(double **A, int N)
{
    int i, j;
    double **res;
    res = makeMatrixOfDoubles(N, N);
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            res[i][j] = A[j][i];
        }
    }
    return res;
}

double **calcT(Eigen *eigenArr, int N, int k)
{
    double **U, rowSum;
    int i, j;
    U = makeMatrixOfDoubles(N, k);
    for (i = 0; i < k; i++)
    {
        for (j = 0; j < N; j++)
        {
            U[j][i] = eigenArr[i].eigenVec[j];
        }
    }
    for (i = 0; i < N; i++)
    {
        rowSum = 0;
        for (j = 0; j < k; j++)
        {
            rowSum += pow(U[i][j], 2);
        }
        rowSum = pow(rowSum, 0.5);
        for (j = 0; j < k; j++)
        {
            if (rowSum==0){
                U[i][j]=0;
            } else {
            U[i][j] /= rowSum;
            }
        }
    }
    return U;
}

/*-----------------NSC------------------*/

/* NSC = Normalized Spectral Clustering */

double **NSC(double **vectorsList, int N, int *d, int k)
{
    double **wam, **Lnorm, **ddm, **V, **T;
    Eigen *eigenArr;
    /*step 1: Generate WAM Matrix*/
    wam = calcWAM(vectorsList, N, *d);
    /*step 2: Compute Lnorm*/
    ddm = calcDDM(wam, N);
    Lnorm = calcLnorm(ddm, wam, N);
    /*step 3: Determine k*/
    V = jacobi(Lnorm, N);
    eigenArr = eigenFill(Lnorm, V, N);
    if (k == 0)
        k = eigenGap(eigenArr, N);
    *d = k;
    /*step 4+5: Generate normalized k eigenVectors Matrix*/
    T = calcT(eigenArr, N, k);
    return T;
}

/*-------------------------END NSC-----------------------*/

/*------------------------------------------Kmeanspp---------------------------------------*/

/* finding the centroids with initial centroids from python:main code*/
void kmeans(int N, int d, int K, int maxiter, double eps, double ***centroidsPointer, double **vectorsList)
{
    int i, itercnt, norm, *cluster_mem;
    double *centroid, **oldcentroids;
    oldcentroids = *centroidsPointer;

    /* int cluster_mem[N]; ~ */
    cluster_mem = calloc(N, sizeof(int));
    if (cluster_mem == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    /* Kmeans code interpeted from pseudo code*/
    itercnt = 0;
    while (1)
    {
        itercnt += 1;
        for (i = 0; i < N; i++)
            cluster_mem[i] = getclosestmeanindex(vectorsList[i], *centroidsPointer, d, K);

        /* double centroids[K][d]; ~ */
        centroid = calloc(K * d, sizeof(double));
        if (centroid == NULL)
        {
            printf("An Error Has Occurred");
            exit(1);
        }
        *centroidsPointer = calloc(K, sizeof(double *));
        if (*centroidsPointer == NULL)
        {
            printf("An Error Has Occurred");
            exit(1);
        }
        for (i = 0; i < K; i++)
        {
            centroidsPointer[0][i] = centroid + i * d;
            calcnewmean(i, d, N, cluster_mem, vectorsList, centroidsPointer[0][i]);
        }

        norm = calcnorm(d, K, eps, *centroidsPointer, oldcentroids);
        free(*oldcentroids);
        free(oldcentroids);
        oldcentroids = *centroidsPointer;
        if (itercnt >= maxiter || norm)
            break;
    }
    free(vectorsList);
    free(cluster_mem);
}

/*Calculate the norm between the old and new centroids. There are K centroids, each has d dimension*/
int calcnorm(int d, int K, double eps, double **centroids, double **oldcentroids)
{
    int i, j;
    double val, *sub;

    /* double sub[d]; ~ */
    sub = calloc(d, sizeof(double));
    if (sub == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    for (i = 0; i < K; i++)
    {
        val = 0.0;
        for (j = 0; j < d; j++)
        {
            sub[j] = centroids[i][j] - oldcentroids[i][j];
            val += sub[j] * sub[j];
        }
        if (sqrt(val) >= eps)
        {
            free(sub);
            return 0;
        }
    }
    free(sub);
    return 1;
}

/* Calculate the centroids after the clustering*/
void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *centroidI)
{
    int k, j, groupsize;
    double *sumvec;
    groupsize = 0;

    /* double sumvec[d] = {0}; ~ */
    sumvec = calloc(d, sizeof(double));
    if (sumvec == NULL)
    {
        printf("An Error Has Occurred");
        exit(1);
    }

    for (k = 0; k < N; k++)
    {
        if (cluster_mem[k] == i)
        {
            for (j = 0; j < d; j++)
            {
                sumvec[j] += vectorsList[k][j];
            }
            groupsize++;
        }
    }
    for (j = 0; j < d; j++)
    {
        centroidI[j] = sumvec[j] / groupsize;
    }
    free(sumvec);
}

/* returns the index of which centroid is the closetst to the vector*/
int getclosestmeanindex(double *vector, double **centroids, int d, int K)
{
    int i, j, minindex;
    double minval;
    minval = 0;
    minindex = 0;
    for (i = 0; i < K; i++)
    {
        double val, *sub;

        /* double sub[d]; ~ */
        sub = calloc(d, sizeof(double));
        if (sub == NULL)
        {
            printf("An Error Has Occurred");
            exit(1);
        }

        val = 0;
        for (j = 0; j < d; j++)
        {
            sub[j] = vector[j] - centroids[i][j];
            val += sub[j] * sub[j];
        }
        if (i == 0)
            minval = val;
        if (val < minval)
        {
            minindex = i;
            minval = val;
        }
        free(sub);
    }
    return minindex;
}

/*-------------------END Kmeanspp-------------------*/

/*-------------------General Functions-------------------*/

int isSymmetric(double **mat, int N)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = i + 1; j < N; j++)
        {
            if (mat[i][j] != mat[j][i])
                return 0;
        }
    }
    return 1;
}

void getNAndD(FILE *ifp, int *NAnddp)
{
    int N, d;
    char c;
    N = 0;
    d = 0;
    while ((c = fgetc(ifp)) != EOF)
    {
        if (c == '\n')
            N++;
        if (c == ',')
            d++;
    }
    NAnddp[0] = N;
    NAnddp[1] = d / N + 1;
}

void makeVeclist(FILE *ifp, double **vectorsList, int N, int d)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (fscanf(ifp, "%lf,", &vectorsList[i][j]) != 1)
            {
                printf("Invalid Input!");
                exit(1);
            }
        }
    }
}

double isNotZero(double dElem)
{
    if (dElem < 0 && dElem > -0.00005)
    {
        return 0;
    }
    return dElem;
}

void printMat(double **mat, int N, int d)
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < d; j++)
        {
            if (j != d - 1)
                printf("%.4f,", mat[i][j]);
            else
                printf("%.4f", mat[i][j]);
        }
        printf("\n");
    }
}

int implementGoalC(int N, int d, char *goal, double **vectorsList)
{
    int i;
    double **wam, **ddg, **lnorm, **V;
    if (strcmp(goal, "wam") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        printMat(wam, N, N);
        free(*wam);
        free(wam);
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        ddg = calcDDM(wam, N);
        printMat(ddg, N, N);
        free(*wam);
        free(wam);
        free(*ddg);
        free(ddg);
    }
    else if (strcmp(goal, "lnorm") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        ddg = calcDDM(wam, N);
        lnorm = calcLnorm(ddg, wam, N);
        printMat(lnorm, N, N);
        free(*wam);
        free(wam);
        free(*ddg);
        free(ddg);
        free(*lnorm);
        free(lnorm);
    }
    else if (strcmp(goal, "jacobi") == 0)
    {
        if (N != d || !isSymmetric(vectorsList, N))
        {
            printf("Invalid Input!");
            return 1;
        }
        V = jacobi(vectorsList, N);
        /*Prints the EigenValues*/
        for (i = 0; i < N; i++)
        {
            if (i != N - 1)
                printf("%.4f,", isNotZero(vectorsList[i][i]));
            else
                printf("%.4f", isNotZero(vectorsList[i][i]));
        }
        printf("\n");
        /*Prints the EigenVectors*/
        printMat(V, N, N);
        free(*V);
        free(V);
    }
    else
    {
        printf("Invalid Input!");
        return 1;
    }
    return 0;
}
double** implementGoalPy(int N, int d, char *goal, double **vectorsList)
{
    int i,j;
    double **wam, **ddg, **lnorm, **V, **jac;
    if (strcmp(goal, "wam") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        return wam;
    }
    else if (strcmp(goal, "ddg") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        ddg = calcDDM(wam, N);
        free(*wam);
        free(wam);
        return ddg;
    }
    else if (strcmp(goal, "lnorm") == 0)
    {
        wam = calcWAM(vectorsList, N, d);
        ddg = calcDDM(wam, N);
        lnorm = calcLnorm(ddg, wam, N);
        free(*wam);
        free(wam);
        free(*ddg);
        free(ddg);
        return lnorm;
    }
    else if (strcmp(goal, "jacobi") == 0)
    {
        if (N != d || !isSymmetric(vectorsList, N))
        {
            printf("Invalid Input!");
            exit(1);
        }
        V = jacobi(vectorsList, N);
        jac=makeMatrixOfDoubles(N+1,N);
        /*jac is matrix that the first row is the  eigen values and the other are the eigenvectors matrix*/
        for (i = 0; i < N; i++) jac[0][i]=vectorsList[i][i];
        for (i=1;i<N+1;i++){
            for (j=0;j<N;j++) jac[i][j]=V[i-1][j];
        }
        free(*V);
        free(V);
        return jac;
    }
    else
    {
        printf("Invalid Input!");
        exit(1);
    }
}


/*-------------------END General Functions-------------------*/

/*-----------------Main Code------------------*/
int main(int argc, char *const argv[]){
    char *goal, *filename;
    double **vectorsList;
    int N, d, *NAndd, exitCode;
    FILE *ifp;
    if (argc != 3)
    {
        printf("Invalid Input!");
        exit(1);
    }
    goal = argv[1];
    filename = argv[2];
    ifp = fopen(filename, "r");
    if (ifp == NULL)
    {
        printf("Invalid Input!");
        exit(1);
    }
    NAndd = calloc(2, sizeof(int));
    assert(NAndd != NULL);
    getNAndD(ifp, NAndd);
    fseek(ifp, 0, SEEK_SET);
    N = NAndd[0];
    d = NAndd[1];
    free(NAndd);
    vectorsList = makeMatrixOfDoubles(N, d);
    makeVeclist(ifp, vectorsList, N, d);
    exitCode = implementGoalC(N, d, goal, vectorsList);
    fclose(ifp);
    free(*vectorsList);
    free(vectorsList);
    exit(exitCode);
    
    /*double **vectorsList,**T_mat,**centroids;
    int  *NAndd,i,d,N;
    FILE *ifp;
    ifp=fopen("input_1.txt","r");
    NAndd = calloc(2, sizeof(int));
    assert(NAndd != NULL);
    getNAndD(ifp, NAndd);
    fseek(ifp, 0, SEEK_SET);
    N = NAndd[0];
    d = NAndd[1];
    free(NAndd);
    vectorsList = makeMatrixOfDoubles(N, d);
    makeVeclist(ifp, vectorsList, N, d);
    T_mat = NSC(vectorsList, N, &d, 0);
    
    centroids=makeMatrixOfDoubles(d,d);
    for (i=0;i<3;i++){
        centroids[0][i]=T_mat[5][i];
    }
    for (i=0;i<3;i++){
        centroids[1][i]=T_mat[8][i];
    }
    for (i=0;i<3;i++){
        centroids[2][i]=T_mat[7][i];
    }
    kmeans(N,d,d,300,0,&centroids,T_mat);
    printMat(centroids,d,d);
    return 0;*/
}




