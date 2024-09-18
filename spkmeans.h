#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define MAXROTATIONS 100
#define EPS 1e-5

typedef struct {
    double c;
    double s;
    int i;
    int j;
} PMat;

typedef struct {
    double eigenVal;
    double *eigenVec;
    int orIndex;
} Eigen;


double ** calcWAM(double **vectorsList, int N, int d);
double calcWeight(double *vec_a, double *vec_b, int d);
double ** calcDDM(double **wam, int N);
double ** calcLnorm(double **ddm, double **wam, int N);
void diagMatMUL(double** mat1,double** mat2,double ***resMatPtr, int N, int whereDiag);
double ** jacobi(double ** A, int N);
void buildP(double ** A, PMat * currP, int N);
void TransformA(double ***APtr,PMat *currP, int N);
void calcNewV(double *** VPtr,PMat *currP, int N);
double ** makeMatrixOfDoubles(int N, int d);
int eigenGap(Eigen *eigenArr, int N);
Eigen * eigenFill(double **Lnorm, double **V, int N);
int comparator(const void *p1, const void *p2);
double ** transpose(double **A, int N);
double ** calcT(Eigen *eigenArr ,int N, int k);
double ** NSC(double ** vectorsList, int N, int* d, int k);
void kmeans(int N, int d, int K, int maxiter, double eps, double ***centroidsPointer, double **vectorsList);
int calcnorm(int d, int K, double eps, double **centroids, double **oldcentroids);
void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *centroidI);
int getclosestmeanindex(double *vector, double **centroids, int d, int K);
int isSymmetric(double **mat, int N);
void getNAndD(FILE *ifp, int *NAnddp);
void makeVeclist(FILE *ifp, double **vectorsList, int N, int d);
double isNotZero(double dElem);
void printMat(double **mat, int N, int d);
int implementGoalC(int N, int d, char *goal, double **vectorsList);
double** implementGoalPy(int N, int d, char *goal, double **vectorsList);
