#define PY_SSIEE_T_CLEAN
#include <Python.h>
#include <math.h>

static void kmeans(int N, int d, int K, int maxiter, double eps, double ***centroidsPointer, double **vectorsList);
static int calcnorm(int d, int K, double eps, double **centroids, double **oldcentroids);
static void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *centroidI);
static int getclosestmeanindex(double *vector, double **centroids, int d, int K);


/*This is the main function to connect C and Python*/
static PyObject *fit(PyObject *self, PyObject *args) {
    double eps, *centroid, **vectorsList, **centroids,*vector;
    int i, j, K, N, d, maxiter;
    PyObject *centroids_from_py, *vectors_from_py, *pyVector, *centroids_to_py;
    if(!PyArg_ParseTuple(args, "iiiidOO", &K, &N, &d, &maxiter, &eps, &vectors_from_py, &centroids_from_py)) {
        return NULL; 
    }
    if (!PyList_Check(centroids_from_py) || !PyList_Check(vectors_from_py)) {
        return NULL;
    }
    /* Allocation: double centroids[K][d];*/
    centroid = calloc(K * d, sizeof(double));
    if (centroid == NULL) {
        return NULL;
    }
    centroids = calloc(K, sizeof(double *));
    if (centroids == NULL) {
        return NULL;
    }
    for (i = 0; i < K; i++) {
        centroids[i] = centroid + i * d;
    }
    
    /* Allocation: double vectorsList[N][d]; ~ */
    vector = calloc(N * d, sizeof(double));
    if (vector == NULL) {
        return NULL;
    }
    vectorsList = calloc(N, sizeof(double *));
    if (vectorsList == NULL) {
        return NULL;
    }
    for (i = 0; i < N; i++) {
        vectorsList[i] = vector + i * d;
    }

    /* Fill vectorsList with the given input vectors from Python */
    for (i = 0; i < N; i++) {
        pyVector = PyList_GetItem(vectors_from_py, i);
        if(!PyList_Check(pyVector)) return NULL;
        for (j = 0; j < d; j++) {
            vectorsList[i][j] = PyFloat_AsDouble(PyList_GetItem(pyVector, j));
            if (PyErr_Occurred() && vectorsList[i][j]  == -1.0){
                if (PyErr_Occurred()) {
                    return NULL;
                }
            }
        } 
    }

    /* Fill intial centroids with the as calculated in Python */
    for (i = 0; i < K; i++) {
        pyVector = PyList_GetItem(centroids_from_py, i);
        if(!PyList_Check(pyVector)) return NULL;
        for (j = 0; j < d; j++) {
            centroids[i][j] = PyFloat_AsDouble(PyList_GetItem(pyVector, j));
            if (PyErr_Occurred() && centroids[i][j]  == -1.0){
                if (PyErr_Occurred()) {
                    return NULL;
                }
            }
        } 
    }
    /*Calculate the desired centroids using Kmeans algorithm*/
    kmeans(N, d, K, maxiter, eps, &centroids, vectorsList);
    free(vector);

    /*Cast the centroids to the python Object and return it*/
    centroids_to_py=PyList_New(K);
    if (centroids_to_py == NULL){
        return NULL;
    }
    for (i = 0; i < K; i++) {
        pyVector = PyList_New(d);
        if (pyVector == NULL) return NULL;
        for (j = 0; j < d; j++) {
            PyList_SetItem(pyVector, j, Py_BuildValue("d", centroids[i][j]));
        }
        PyList_SetItem(centroids_to_py, i, Py_BuildValue("O", pyVector));
    }
    free(centroids);
    return centroids_to_py;
}


/* finding the centroids with initial centroids from python:main code*/
static void kmeans(int N, int d, int K, int maxiter, double eps, double ***centroidsPointer, double **vectorsList) {
    int i, itercnt, norm, *cluster_mem;
    double *centroid, **oldcentroids;
    oldcentroids = *centroidsPointer;

    /* int cluster_mem[N]; ~ */
    cluster_mem = calloc(N, sizeof(int));
    if (cluster_mem == NULL) {
        printf("An Error Has Occurred");
        exit(1);
    }

/* Kmeans code interpeted from pseudo code*/
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

        norm = calcnorm(d, K, eps, *centroidsPointer, oldcentroids);
        free(*oldcentroids); free(oldcentroids);
        oldcentroids = *centroidsPointer;
        if (itercnt >= maxiter || norm) break;
    }
    free(vectorsList); free(cluster_mem);
}

/*Calculate the norm between the old and new centroids. There are K centroids, each has d dimension*/
static int calcnorm(int d, int K, double eps, double **centroids, double **oldcentroids) {
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
        if (sqrt(val) >= eps) {
            free(sub);
            return 0;
        }
    }
    free(sub);
    return 1;
}

/* Calculate the centroids after the clustering*/
static void calcnewmean(int i, int d, int N, int *cluster_mem, double **vectorsList, double *centroidI) {
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

/* returns the index of which centroid is the closetst to the vector*/
static int getclosestmeanindex(double *vector, double **centroids, int d, int K) {
    int i, j, minindex;
    double minval;
    minval = 0;
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

/*C-Python-Functions*/
static PyMethodDef kmeansMethods[] = {
        {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Implements the Kmeans algorithm with C")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        kmeansMethods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}
