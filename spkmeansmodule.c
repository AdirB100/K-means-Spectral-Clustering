#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject *notSpk(PyObject *self, PyObject *args);
static PyObject *getTmatPy(PyObject *self, PyObject *args);
static PyObject *fit(PyObject *self, PyObject *args);

double **T_mat;
int N,d;

/**for everything but spk, according to given guides in pdf, should print the desired value and return 0*/
static PyObject *notSpk(PyObject *self, PyObject *args) {
    PyObject *vectors_from_py, *pyVector,*mat_to_py;
    int i, j, k,len;
    double **vectorsList, **pymat;
    char *goal;
    if (!PyArg_ParseTuple(args, "iiisO", &N, &d, &k, &goal, &vectors_from_py)) {
        return NULL;
    }

    if (!PyList_Check(vectors_from_py)) {
        return NULL;
    }

    vectorsList = makeMatrixOfDoubles(N, d);

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

    pymat= implementGoalPy(N, d, goal, vectorsList);
    free(*vectorsList); free(vectorsList);
    if (pymat==NULL) return NULL;
    if (strcmp(goal,"jacobi")==0) len=N+1;
    else len=N;
    mat_to_py=PyList_New(len);
    if (mat_to_py == NULL){
        return NULL;
    }
    for (i = 0; i < len; i++) {
        pyVector = PyList_New(N);
        if (pyVector == NULL) return NULL;
        for (j = 0; j < N; j++) {
            PyList_SetItem(pyVector, j, Py_BuildValue("d", pymat[i][j]));
        }
        PyList_SetItem(mat_to_py, i, Py_BuildValue("O", pyVector));
    }
    return mat_to_py;
}



/**initializing all needed matrices for spk algorithm*/
static PyObject *getTmatPy(PyObject *self, PyObject *args){
    PyObject *vectors_from_py, *pyVector, *T_mat_to_py;
    int i, j, k;
    double **vectorsList;

    
    if (!PyArg_ParseTuple(args, "iiiO", &N, &d, &k, &vectors_from_py)) {
        return NULL;
    }
    if (!PyList_Check(vectors_from_py)) {
        return NULL;
    }
    
    vectorsList = makeMatrixOfDoubles(N, d);

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

    
    T_mat = NSC(vectorsList, N, &d, k);

    /*Cast the centroids to the python Object and return it*/
    T_mat_to_py=PyList_New(N);
    if (T_mat_to_py == NULL){
        return NULL;
    }
    for (i = 0; i < N; i++) {
        pyVector = PyList_New(d);
        if (pyVector == NULL) return NULL;
        for (j = 0; j < d; j++) {
            PyList_SetItem(pyVector, j, Py_BuildValue("d", T_mat[i][j]));
        }
        PyList_SetItem(T_mat_to_py, i, Py_BuildValue("O", pyVector));
    }

    return T_mat_to_py;
}

static PyObject *fit(PyObject *self, PyObject *args) {
    double **centroids;
    int i, j, K;
    K=d;
    PyObject *centroids_from_py, *pyVector,*centroids_to_py;
    if(!PyArg_ParseTuple(args, "O", &centroids_from_py)) {
        return NULL; 
    }
    if (!PyList_Check(centroids_from_py)) {
        return NULL;
    }

    centroids = makeMatrixOfDoubles(K, d);
    
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
    kmeans(N, d, K, 300, 0, &centroids, T_mat);

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
    return centroids_to_py;
}

/*module:*/
static PyMethodDef kmeansMethods[] = {
        {"notSpk", (PyCFunction) notSpk, METH_VARARGS, PyDoc_STR("Does the goals which different from spk")},
        {"getTmatPy", (PyCFunction) getTmatPy, METH_VARARGS, PyDoc_STR("Finds T matrix")},
        {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Finds the centroids with kmeans algorithm")},
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