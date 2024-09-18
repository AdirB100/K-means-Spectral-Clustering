#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

void kmeans(int N, int d, int K, FILE *ifp, int maxiter, double ***centroidsPointer);
int calcnorm(int d, int K, double **newcentroids, double **centroids);
void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *newcentroid);
int getclosestmeanindex(double *vector, double **centroids, int d, int K);
void makeVeclist(FILE *ifp, double **vectorsList, int N, int d);
void getNAndD(FILE *ifp, int *NAnddp);

int main(int argc, char *argv[]) {
    int N, d, K, i, j, maxiter, *NAndd;
    unsigned m;
    double *centroid, **centroids;
    char *ifp_path, *ofp_path;
    FILE *ifp, *ofp;
    if (argc != 5 && argc != 4) {
        printf("Invalid Input!");
        exit(1);
    }
    for (m = 0; m < strlen(argv[1]); m++) {
        if (!isdigit(argv[1][m])) {
            printf("Invalid Input!");
            exit(1);
        }
    }
    K = atoi(argv[1]);
    if (K <= 1) {
        printf("Invalid Input!");
        exit(1);
    }
    if (argc < 5) {
        maxiter = 200;
        ifp_path = argv[2];
        ofp_path = argv[3];
    } else {
        for (m = 0; m < strlen(argv[2]); m++) {
            if (!isdigit(argv[2][m])) {
                printf("Invalid Input!");
                exit(1);
            }
        }
        maxiter = atoi(argv[2]);
        if (maxiter == 0) {
            printf("Invalid Input!");
            exit(1);
        }
        ifp_path = argv[3];
        ofp_path = argv[4];
    }
    ifp = fopen(ifp_path, "r");
    if (ifp == NULL) {
        printf("Invalid Input!");
        exit(1);
    }

    /* int NAndd[2]; ~ */
    NAndd = calloc(2, sizeof(int));
    if (NAndd == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }

    getNAndD(ifp, NAndd);
    fseek(ifp, 0, SEEK_SET);
    N = NAndd[0]; d = NAndd[1];
    free(NAndd);
    if (K >= N) {
        printf("Invalid Input!");
        exit(1);
    }
    
    /* double centroids[K][d]; ~ */
    centroid = calloc(K * d, sizeof(double));
    if (centroid == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    centroids = calloc(K, sizeof(double *));
    if (centroids == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < K; i++) {
        centroids[i] = centroid + i * d;
    }

    kmeans(N, d, K, ifp, maxiter, &centroids);
    fclose(ifp);
    ofp = fopen(ofp_path, "w");
    if (ofp == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    for (i = 0; i < K; i++) {
        for (j = 0; j < d; j++) {
            fprintf(ofp, "%.4f", centroids[i][j]);
            if (j != d - 1) fputc(',', ofp);
        }
        fputc('\n', ofp);
    }
    fclose(ofp);
    free(*centroids); free(centroids);
    return 0;
}

void kmeans(int N, int d, int K, FILE *ifp, int maxiter, double ***centroidsPointer) {
    int i, j, itercnt, norm, *cluster_mem;
    double *vector, *centroid, **vectorsList, **oldcentroids;

    /* double vectorsList[N][d]; ~ */
    vector = calloc(N * d, sizeof(double));
    if (vector == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    vectorsList = calloc(N, sizeof(double *));
    if (vectorsList == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    
    for (i = 0; i < N; i++) vectorsList[i] = vector + i * d;
    oldcentroids = *centroidsPointer;

    /* int cluster_mem[N]; ~ */
    cluster_mem = calloc(N, sizeof(int));
    if (cluster_mem == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }

    makeVeclist(ifp, vectorsList, N, d);
    for (i = 0; i < K; i++) for (j = 0; j < d; j++) centroidsPointer[0][i][j] = vectorsList[i][j];
    itercnt = 0;
    while (1) {
        itercnt += 1;
        for (i = 0; i < N; i++) cluster_mem[i] = getclosestmeanindex(vectorsList[i], *centroidsPointer, d, K);

        /* double centroids[K][d]; ~ */
        centroid = calloc(K * d, sizeof(double));
        if (centroid == NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }
        *centroidsPointer = calloc(K, sizeof(double *));
        if (*centroidsPointer == NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }
        for (i = 0; i < K; i++) {
            centroidsPointer[0][i] = centroid + i * d;
            calcnewmean(i, d, N, cluster_mem, vectorsList, centroidsPointer[0][i]);
        }

        norm = calcnorm(d, K, *centroidsPointer, oldcentroids);
        free(*oldcentroids); free(oldcentroids);
        oldcentroids = *centroidsPointer;
        if (itercnt >= maxiter || norm) break;
    }
    free(vector); free(vectorsList); free(cluster_mem);
}

int calcnorm(int d, int K, double **centroids, double **oldcentroids) {
    int i, j;
    double val, *sub;

    /* double sub[d]; ~ */
    sub = calloc(d, sizeof(double));
    if (sub == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }
    
    for (i = 0; i < K; i++) {
        val = 0.0;
        for (j = 0; j < d; j++) {
            sub[j] = centroids[i][j] - oldcentroids[i][j];
            val += sub[j] * sub[j];
        }
        if (sqrt(val) >= 0.001) {
            free(sub);
            return 0;
        }
    }
    free(sub);
    return 1;
}

void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *centroidI) {
    int k, j, groupsize;
    double *sumvec;
    groupsize = 0;

    /* double sumvec[d] = {0}; ~ */
    sumvec = calloc(d, sizeof(double));
    if (sumvec == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }

    for (k = 0; k < N; k++) {
        if (cluster_mem[k] == i) {
            for (j = 0; j < d; j++) {
                sumvec[j] += vectorsList[k][j];
            }
            groupsize++;
        }
    }
    for (j = 0; j < d; j++) {
        centroidI[j] = sumvec[j] / groupsize;
    }
    free(sumvec);
}

int getclosestmeanindex(double *vector, double **centroids, int d, int K) {
    int i, j, minindex;
    double minval;

    minindex = 0;
    for (i = 0; i < K; i++) {
        double val, *sub;

        /* double sub[d]; ~ */
        sub = calloc(d, sizeof(double));
        if (sub == NULL) {
            printf("An Error Has Occurred");
            exit(1);
        }

        val = 0;
        for (j = 0; j < d; j++) {
            sub[j] = vector[j] - centroids[i][j];
            val += sub[j] * sub[j];
        }
        if (i == 0) minval = val;
        if (val < minval) {
            minindex = i;
            minval = val;
        }
        free(sub);
    }
    return minindex;
}

void makeVeclist(FILE *ifp, double **vectorsList, int N, int d) {
    int i, j;
    for (i = 0; i < N; i++) for (j = 0; j < d; j++) fscanf(ifp, "%lf,", &vectorsList[i][j]);
}

void getNAndD(FILE *ifp, int *NAnddp) {
    int N, d;
    char c;
    N = 0; d = 0;
    while ((c = fgetc(ifp)) != EOF) {
        if (c == '\n') N++;
        if (c == ',') d++;
    }
    NAnddp[0] = N;
    NAnddp[1] = d / N + 1;
}
