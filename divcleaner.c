#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#define PyArray_FSimpleNewFromData(nd, dims, typenum, data) \
        PyArray_New(&PyArray_Type, nd, dims, typenum, NULL, \
                    data, 0, NPY_ARRAY_FARRAY, NULL)

#include <Python.h>
#include <numpy/arrayobject.h>
#include "bicgstab.h"

/* Docstrings */

static char module_docstring[] =
  "This modules provides an interface for solving the Poisson equation.";
static char divcleaner_docstring[] =
  "Solve the Poisson equation";
static char phiToBbc_docstring[] =
  "Convert field potential to boundary-centred magnetic field, assuming \
	  periodic boundary conditions in all directions.";
static char fixTopBottomBoundary_docstring[] =
  "Fix the values of top and bottom boundaries of vertical magnetic field \
	  assuming div(B) = 0.";
static char div_docstring[] =
  "Compute divergence of vector field.";
static char killdiv_docstring[] =
  "Remove the irrotational part of a given field.";

/* Module function prototypes */

static PyObject* divcleaner_poisson(PyObject *self, PyObject *args, \
		PyObject *kwargs);

static PyObject* divcleaner_phiToBbc(PyObject *self, PyObject *args, \
		PyObject *kwargs);

static PyObject* divcleaner_fixTopBottomBoundary(PyObject *self, \
		PyObject *args, PyObject *kwargs);

static PyObject* divcleaner_div(PyObject *self, PyObject *args, \
		PyObject *kwargs);

static PyObject* divcleaner_killdiv(PyObject *self, PyObject *args, \
		PyObject *kwargs);

/* Module methods declaration */

static PyMethodDef module_methods[] = {
	{"poisson", (PyCFunction)divcleaner_poisson, \
		METH_VARARGS | METH_KEYWORDS, divcleaner_docstring},
	{"phiToBbc", (PyCFunction)divcleaner_phiToBbc, \
		METH_VARARGS | METH_KEYWORDS, phiToBbc_docstring},
	{"fixTopBottomBoundary", \
		(PyCFunction)divcleaner_fixTopBottomBoundary, \
		METH_VARARGS | METH_KEYWORDS, fixTopBottomBoundary_docstring},
	{"div", (PyCFunction)divcleaner_div, \
		METH_VARARGS | METH_KEYWORDS, div_docstring},
	{"killdiv", (PyCFunction)divcleaner_killdiv, \
		METH_VARARGS | METH_KEYWORDS, killdiv_docstring},
	{NULL, NULL, 0, NULL}
};

static char* divcleaner_kwlist[] = {
	"x0", "bc", "thr", "maxiter",
	"b", "boundary", "grid_spacing",
	NULL
};

static char* phiToBbc_kwlist[] = {
	"phi", "grid_spacing",
	NULL
};

static char* fixTopBottomBoundary_kwlist[] = {
	"Bb1", "Bb2", "Bb3",
	"grid_spacing",
	NULL
};

static char* div_kwlist[] = {
	"Bb1", "Bb2", "Bb3",
	"grid_spacing",
	NULL
};

static char* killdiv_kwlist[] = {
	"Bb1", "Bb2", "Bb3",
	"grid_spacing", "bc", "boundary",
	"x0", "thr", "maxiter",
	NULL
};

/* Module initialization */

PyMODINIT_FUNC initdivcleaner(void) {
	PyObject *m = Py_InitModule3("divcleaner", module_methods, \
			module_docstring);
	if (m == NULL) {
		/*
		 PyErr_SetString(PyExc_ValueError, \
				 "initdivcleaner: Failed loading module.");
		*/
		return;
	}

	/* Load `numpy` functionality */
	import_array();
}

static PyObject* divcleaner_poisson(PyObject *self, PyObject *args, \
		PyObject *kwargs) {
	/* Input arguments */
	PyObject *x0_obj = Py_BuildValue("s", NULL);
	int bX, bY, bZ;
	double thr;
	int maxiter;

	/* Keyword arguments */
	PyObject *b_obj  = Py_BuildValue("s", NULL), \
		 *d_obj  = Py_BuildValue("s", NULL);
	double hx=1., hy=1., hz=1.;

	/* Output variables */
	PyObject *x_arr  = Py_BuildValue("s", NULL);

	/* Local variables */
	int ndim, a1, a2, a3, n, N, nghost=1;
	void* stream = NULL;

	double *b  = NULL, \
	       *x0 = NULL, \
	       *S  = NULL, \
	       *x  = NULL, \
	       *g  = NULL, \
	       *x_out = NULL;

	PyArrayObject *b_arr  = NULL, \
	              *x0_arr = NULL, \
		      *d_arr  = NULL;

	npy_intp *dims = NULL, \
	         *dims_d = NULL, \
	         *dims_out = NULL;
	
	Lambda L = (Lambda)&stencil3DL;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"O(iii)di|OO(ddd)", divcleaner_kwlist, \
				&x0_obj, &bX, &bY, &bZ, &thr, &maxiter, \
				&b_obj, &d_obj, &hx, &hy, &hz)) {
		return NULL;
	}

	/* Interprete input objects as numpy arrays */
	x0_arr = (PyArrayObject*)PyArray_FROM_OTF(x0_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	if (b_obj != Py_None)
		b_arr = (PyArrayObject*)PyArray_FROM_OTF(b_obj, \
				NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	if (d_obj == Py_None && \
			(bX == BC_DIRICHLET || bY == BC_DIRICHLET \
			 || bZ == BC_DIRICHLET)) {
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_poisson: With Dirichlet boundary \
				conditions the keyword argument `boundary` \
				must be provided.");
		return NULL;
	}
	if (d_obj != Py_None) {
		d_arr = (PyArrayObject*)PyArray_FROM_OTF(d_obj, \
				NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	}

	/* Throw an exception if that didn't work */
	if (x0_arr == NULL || (b_arr == NULL && b_obj != Py_None)) {
		Py_XDECREF(x0_arr);
		Py_XDECREF(b_arr);
		Py_XDECREF(d_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_poisson: Failed parsing `x0` \
				and/or `b`.");
		return NULL;
	}

	/* Get dimensions of the computational box */
	ndim = PyArray_NDIM(x0_arr);
	if (ndim != 3) {
		Py_XDECREF(x0_arr);
		Py_XDECREF(b_arr);
		Py_XDECREF(d_arr);
		PyErr_SetString(PyExc_ValueError, "divcleaner_poisson: \
				`x0` does not have a valid shape.");
		return NULL;
	}
	dims = PyArray_DIMS(x0_arr);
	a1 = (int)dims[0]+2*nghost;
	a2 = (int)dims[1]+2*nghost;
	a3 = (int)dims[2]+2*nghost;
	dims_out = PyMem_Malloc(3*sizeof(npy_intp));
	dims_out[0] = a1;
	dims_out[1] = a2;
	dims_out[2] = a3;
	n = dims[0]*dims[1]*dims[2];
	N = a1*a2*a3;

	/* Throw exception if b or boundary do not have a valid shape */
	if (b_obj != Py_None) {
		ndim = PyArray_NDIM(b_arr);
		if (ndim != 3) {
			PyErr_SetString(PyExc_ValueError, \
					"divcleaner_poisson: `b` argument \
					does not have a valid shape.");
			Py_XDECREF(x0_arr);
			Py_XDECREF(b_arr);
			Py_XDECREF(d_arr);
			return NULL;
		}
		dims = PyArray_DIMS(b_arr);
		if ((int)dims[0]+2*nghost != a1 || (int)dims[1]+2*nghost != a2
				|| (int)dims[2]+2*nghost != a3) {
			PyErr_SetString(PyExc_ValueError, \
					"divcleaner_poisson: `b` argument \
					does not have a valid shape.");
			Py_XDECREF(x0_arr);
			Py_XDECREF(b_arr);
			Py_XDECREF(d_arr);
			return NULL;
		}
	
	}
	if (d_obj != Py_None) {
		ndim = PyArray_NDIM(d_arr);
		if (ndim != 3) {
			PyErr_SetString(PyExc_ValueError, \
					"divcleaner_poisson: `b` argument \
					does not have a valid shape.");
			Py_XDECREF(x0_arr);
			Py_XDECREF(b_arr);
			Py_XDECREF(d_arr);
			return NULL;
		}
		dims_d = PyArray_DIMS(d_arr);
		if ((int)dims_d[0] != a1 || (int)dims_d[1] != a2
				|| (int)dims_d[2] != a3) {
			PyErr_SetString(PyExc_ValueError, \
					"divcleaner_poisson: `boundary` \
					argument does not have a valid \
					shape.");
			Py_XDECREF(x0_arr);
			Py_XDECREF(b_arr);
			Py_XDECREF(d_arr);
			return NULL;
		}
	}

	/* Create the stencil */
	S = stencil2OLaplacian3D(hx, hy, hz);

	/* Create the stream for the laplacian lambda operator */
	if (d_obj != Py_None)
		g = saveGhostCells3D(PyArray_DATA(d_arr), a1, a2, a3, nghost);
	stream = streamForLaplacian3DL(a1, a2, a3, nghost, \
			bX, bY, bZ, S, g);
	free(g);

	/* Allocate C-arrays */
	x0    = arrmalloc(N);
	b     = arrmalloc(N);
	x     = arrmalloc(N);
	x_out = PyMem_Malloc(n*sizeof(double));

	/* Set b */
	if (b_obj != Py_None) {
		nullifyGhostFrame3D(b, a1, a2, a3, nghost);
		addGhostFrame3D(b, PyArray_DATA(b_arr), a1, a2, a3, nghost);
	}
	else nullify(b, N);

	/* Update b taking into account constant boundary conditions */
	nonLinearPart(x, L, a1, a2, a3, nghost, stream);
	add(b, b, -1., x, N);

	/* Update stream */
	if (bX == BC_DIRICHLET) bX = BC_ZERO;
	if (bY == BC_DIRICHLET) bY = BC_ZERO;
	if (bZ == BC_DIRICHLET) bZ = BC_ZERO;
	free(stream);
	stream = streamForLaplacian3DL(a1, a2, a3, nghost, \
			bX, bY, bZ, S, NULL);

	/* Set x0 */
	nullifyGhostFrame3D(x0, a1, a2, a3, nghost);
	addGhostFrame3D(x0, PyArray_DATA(x0_arr), a1, a2, a3, nghost);

	/* Solve the Poisson equation */
	solveL(x, L, stream, b, x0, thr, maxiter, N);

	/* Remove ghost cells */
	removeGhostFrame3D(x_out, x, a1, a2, a3, nghost);
	x_arr = PyArray_FSimpleNewFromData(ndim, dims, NPY_DOUBLE, x_out);

	/* Free memory */	
	Py_XDECREF(x0_arr);
	Py_XDECREF(b_arr);
	Py_XDECREF(d_arr);
	free(b);
	free(x0);
	free(stream);
	free(S);

	return x_arr;
}

static PyObject* divcleaner_phiToBbc(PyObject *self, PyObject *args, \
		PyObject *kwargs) {
	/* Input arguments */
	PyObject *phi_obj = Py_BuildValue("s", NULL);

	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Output variables */
	PyObject *Bb1_arr = Py_BuildValue("s", NULL), \
	         *Bb2_arr = Py_BuildValue("s", NULL), \
	         *Bb3_arr = Py_BuildValue("s", NULL);

	/* Local variables */
	int a1, a2, a3, n1, n2, n3, i, j, k;
	double *phi, *Bb1, *Bb2, *Bb3;
	PyArrayObject *phi_arr  = NULL;
	npy_intp *dims = NULL, *dims1 = NULL, *dims2 = NULL, *dims3 = NULL;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"O|(ddd)", phiToBbc_kwlist, \
				&phi_obj, &hx, &hy, &hz)) {
		return NULL;
	}

	/* Interprete input objects as numpy arrays */
	phi_arr = (PyArrayObject*)PyArray_FROM_OTF(phi_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

	/* Throw an exception if that didn't work */
	if (phi_arr == NULL) {
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_phiToBbc: Failed parsing `phi`.");
		Py_XDECREF(phi_arr);
		return NULL;
	}

	/* Get dimensions of the computational box */
	if (PyArray_NDIM(phi_arr) != 3) {
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_phiToBbc: `phi` does not have \
				a valid shape.");
		Py_XDECREF(phi_arr);
		return NULL;
	}
	dims = PyArray_DIMS(phi_arr);
	a1 = (int)dims[0];
	a2 = (int)dims[1];
	a3 = (int)dims[2];
	n1 = (a1+1)*a2*a3;
	n2 = a1*(a2+1)*a3;
	n3 = a1*a2*(a3+1);
	dims1 = PyMem_Malloc(3*sizeof(npy_intp));
	dims2 = PyMem_Malloc(3*sizeof(npy_intp));
	dims3 = PyMem_Malloc(3*sizeof(npy_intp));
	dims1[0] = a1+1;
	dims1[1] = a2;
	dims1[2] = a3;
	dims2[0] = a1;
	dims2[1] = a2+1;
	dims2[2] = a3;
	dims3[0] = a1;
	dims3[1] = a2;
	dims3[2] = a3+1;

	/* Allocate memory for output */
	Bb1 = PyMem_Malloc(n1*sizeof(double));
	Bb2 = PyMem_Malloc(n2*sizeof(double));
	Bb3 = PyMem_Malloc(n3*sizeof(double));

	/* Compute Bb */
	phi = PyArray_DATA(phi_arr);
	
	for (k=a3-1; k>=0; k--) {
	  for (j=a2-1; j>=0; j--) {
	    for (i=a1-1; i>=1; i--) {
	      Bb1[(k*a2+j)*(a1+1)+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[(k*a2+j)*a1+i-1])/hx;
	    }
	    Bb1[(k*a2+j)*(a1+1)] = \
	         (phi[(k*a2+j)*a1] - phi[(k*a2+j)*a1+a1-1])/hx;
	    Bb1[(k*a2+j)*(a1+1)+a1] = Bb1[(k*a2+j)*(a1+1)];
	  }
	}
	for (k=a3-1; k>=0; k--) {
	  for (i=a1-1; i>=0; i--) {
	    for (j=a2-1; j>=1; j--) {
	      Bb2[(k*(a2+1)+j)*a1+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[(k*a2+j-1)*a1+i])/hy;
	    }
	    Bb2[(k*(a2+1))*a1+i] = \
	         (phi[(k*a2)*a1+i] - phi[(k*a2+a2-1)*a1+i])/hy;
	    Bb2[(k*(a2+1)+a2)*a1+i] = Bb2[(k*(a2+1))*a1+i];
	  }
	}
	for (j=a2-1; j>=0; j--) {
	  for (i=a1-1; i>=0; i--) {
	    for (k=a3-1; k>=1; k--) {
	      Bb3[(k*a2+j)*a1+i] = \
	           (phi[(k*a2+j)*a1+i] - phi[((k-1)*a2+j)*a1+i])/hz;
	    }
	    Bb3[j*a1+i] = \
	         (phi[j*a1+i] - phi[((a3-1)*a2+j)*a1+i])/hz;
	    Bb3[(a3*a2+j)*a1+i] = Bb3[(k*a2+j)*a1+i];
	  }
	}

	Bb1_arr = PyArray_FSimpleNewFromData(3, dims1, NPY_DOUBLE, Bb1);
	Bb2_arr = PyArray_FSimpleNewFromData(3, dims2, NPY_DOUBLE, Bb2);
	Bb3_arr = PyArray_FSimpleNewFromData(3, dims3, NPY_DOUBLE, Bb3);

	return Py_BuildValue("(OOO)", Bb1_arr, Bb2_arr, Bb3_arr);
}

static PyObject* divcleaner_fixTopBottomBoundary(PyObject *self, \
		PyObject *args, PyObject *kwargs) {
	/* Input arguments */
	PyObject *Bb1_obj = Py_BuildValue("s", NULL), \
	         *Bb2_obj = Py_BuildValue("s", NULL), \
	         *Bb3_obj = Py_BuildValue("s", NULL);

	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Output variables */
	PyObject *Bb3_out_arr = NULL;

	/* Local variables */
	int a1, a2, a3, i, j;
	double *Bb1, *Bb2, *Bb3, *Bb3_out;
	PyArrayObject *Bb1_arr = NULL, \
	              *Bb2_arr = NULL, \
	              *Bb3_arr = NULL;
	npy_intp *dims = NULL, *dims_tmp = NULL;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"OOO|(ddd)", fixTopBottomBoundary_kwlist, \
				&Bb1_obj, &Bb2_obj, &Bb3_obj, \
				&hx, &hy, &hz)) {
		return NULL;
	}

	/* Interprete input objects as numpy arrays */
	Bb1_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb1_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb2_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb2_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb3_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb3_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

	/* Throw an exception if that didn't work */
	if (Bb1_arr == NULL || Bb2_arr == NULL || Bb3_arr == NULL) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_fixTopBottomBoundary: Failed \
				parsing `Bb1` and/or `Bb2` and/or `Bb3`.");
		return NULL;
	}

	/* Get dimensions of the computational box */
	if (PyArray_NDIM(Bb1_arr) != 3 || PyArray_NDIM(Bb2_arr) != 3 || \
			PyArray_NDIM(Bb3_arr) != 3) {
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_fixTopBottomBoundary: `Bb1` \
				and/or `Bb2` and/or `Bb3` do not have a \
				valid shape.");
		return NULL;
	}
	dims = PyArray_DIMS(Bb3_arr);

	a1 = dims[0];
	a2 = dims[1];
	a3 = dims[2]-1;

	/* Check if shapes are consistent */
	dims_tmp = PyArray_DIMS(Bb1_arr);
	if (dims_tmp[0] != a1+1 || dims_tmp[1] != a2 || dims_tmp[2] != a3) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_fixTopBottomBoundary: `Bb1` \
				and `Bb3` do not have consistent shapes.");
		return NULL;
	}
	dims_tmp = PyArray_DIMS(Bb2_arr);
	if (dims_tmp[0] != a1 || dims_tmp[1] != a2+1 || dims_tmp[2] != a3) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_fixTopBottomBoundary: `Bb2` \
				and `Bb3` do not have consistent shapes.");
		return NULL;
	}

	/* Allocate memory for output */
	Bb3_out = PyMem_Malloc(a1*a2*(a3+1)*sizeof(double));

	Bb1 = PyArray_DATA(Bb1_arr);
	Bb2 = PyArray_DATA(Bb2_arr);
	Bb3 = PyArray_DATA(Bb3_arr);

	/* Fix top and bottom boundaries */
	arrcpy(Bb3_out, Bb3, a1*a2*(a3+1));
	for (j=a2-1; j>=0; j--) {
	  for (i=a1-1; i>=0; i--) {
	    Bb3_out[j*a1+i] =         - hy*hz*Bb1[j*(a1+1)+i] \
	                              + hy*hz*Bb1[j*(a1+1)+i+1] \
	                              - hz*hx*Bb2[j*a1+i] \
	                              + hz*hx*Bb2[(j+1)*a1+i] \
	                              + hx*hy*Bb3[(a2+j)*a1+i];
	    Bb3_out[j*a1+i] /= hx*hy;
	    Bb3_out[(a3*a2+j)*a1+i] = + hy*hz*Bb1[((a3-1)*a2+j)*(a1+1)+i] \
	                              - hy*hz*Bb1[((a3-1)*a2+j)*(a1+1)+i+1] \
	                              + hz*hx*Bb2[((a3-1)*(a2+1)+j)*a1+i] \
	                              - hz*hx*Bb2[((a3-1)*(a2+1)+j+1)*a1+i] \
	                              + hx*hy*Bb3[((a3-1)*a2+j)*a1+i];
	    Bb3_out[(a3*a2+j)*a1+i] /= hx*hy;
	  }
	}

	/* Free memory */
	Py_XDECREF(Bb1_arr);
	Py_XDECREF(Bb2_arr);
	Py_XDECREF(Bb3_arr);

	/* Construct output */
	Bb3_out_arr = PyArray_FSimpleNewFromData(3, dims, NPY_DOUBLE, Bb3_out);

	return Bb3_out_arr;
}

static PyObject* divcleaner_div(PyObject *self, PyObject *args, \
		PyObject *kwargs) {
	/* Input arguments */
	PyObject *Bb1_obj = Py_BuildValue("s", NULL), \
	         *Bb2_obj = Py_BuildValue("s", NULL), \
	         *Bb3_obj = Py_BuildValue("s", NULL);

	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;

	/* Output variables */
	PyObject *div_arr = NULL;

	/* Local variables */
	int a1, a2, a3, i, j, k;
	double *Bb1, *Bb2, *Bb3, *div;
	PyArrayObject *Bb1_arr = NULL, \
	              *Bb2_arr = NULL, \
	              *Bb3_arr = NULL;
	npy_intp *dims = NULL, *dims_tmp = NULL, *dims_out = NULL;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"OOO|(ddd)", div_kwlist, \
				&Bb1_obj, &Bb2_obj, &Bb3_obj, \
				&hx, &hy, &hz)) {
		return NULL;
	}

	/* Interprete input objects as `numpy` arrays */
	Bb1_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb1_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb2_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb2_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb3_arr = (PyArrayObject*)PyArray_FROM_OTF(Bb3_obj, \
			NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

	/* Throw an exception if that didn't work */
	if (Bb1_arr == NULL || Bb2_arr == NULL || Bb3_arr == NULL) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_div: Failed parsing `Bb1` \
				and/or `Bb2` and/or `Bb3`.");
		return NULL;
	}

	/* Get dimensions of the computational box */
	if (PyArray_NDIM(Bb1_arr) != 3 || PyArray_NDIM(Bb2_arr) != 3 || \
			PyArray_NDIM(Bb3_arr) != 3) {
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_div: `Bb1` and/or `Bb2` and/or \
				`Bb3` do not have a valid shape.");
		return NULL;
	}
	dims = PyArray_DIMS(Bb3_arr);
	a1 = dims[0];
	a2 = dims[1];
	a3 = dims[2]-1;

	/* Check if shapes are consistent */
	dims_tmp = PyArray_DIMS(Bb1_arr);
	if (dims_tmp[0] != a1+1 || dims_tmp[1] != a2 || dims_tmp[2] != a3) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_div: `Bb1` and `Bb3` do not \
				have consistent shapes.");
		return NULL;
	}
	dims_tmp = PyArray_DIMS(Bb2_arr);
	if (dims_tmp[0] != a1 || dims_tmp[1] != a2+1 || dims_tmp[2] != a3) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_div: `Bb2` and `Bb3` do not \
				have consistent shapes.");
		return NULL;
	}

	/* Allocate memory for output */
	div = PyMem_Malloc(a1*a2*a3*sizeof(double));

	Bb1 = PyArray_DATA(Bb1_arr);
	Bb2 = PyArray_DATA(Bb2_arr);
	Bb3 = PyArray_DATA(Bb3_arr);

	/* Compute divergence */
	for (k=a3-1; k>=0; k--) {
	  for (j=a2-1; j>=0; j--) {
	    for (i=a1-1; i>=0; i--) {
	      div[(k*a2+j)*a1+i] = - Bb1[(k*a2+j)*(a1+1)+i]/hx \
	                           + Bb1[(k*a2+j)*(a1+1)+i+1]/hx \
	                           - Bb2[(k*(a2+1)+j)*a1+i]/hy \
	                           + Bb2[(k*(a2+1)+j+1)*a1+i]/hy \
	                           - Bb3[(k*a2+j)*a1+i]/hz \
	                           + Bb3[((k+1)*a2+j)*a1+i]/hz;
	    }
	  }
	}

	/* Free memory */
	Py_XDECREF(Bb1_arr);
	Py_XDECREF(Bb2_arr);
	Py_XDECREF(Bb3_arr);

	/* Construct output */
	dims_out = PyMem_Malloc(3*sizeof(npy_intp));
	dims_out[0] = a1;
	dims_out[1] = a2;
	dims_out[2] = a3;
	div_arr = PyArray_FSimpleNewFromData(3, dims_out, NPY_DOUBLE, div);

	return div_arr;
}

static PyObject* divcleaner_killdiv(PyObject *self, PyObject *args, \
		PyObject *kwargs) {
	/* Input arguments */
	PyObject *Bb1_obj, *Bb2_obj, *Bb3_obj;

	/* Keyword arguments */
	double hx=1., hy=1., hz=1.;
	int bX=BC_PERIODIC, bY=BC_PERIODIC, bZ=BC_ZERO;
	PyObject *boundary=Py_BuildValue("s", NULL), \
	         *x0_obj=Py_BuildValue("s", NULL);
	double thr=1.e-5;
	int maxiter=20;

	/* Output variables */
	PyObject *Bb1_out, *Bb2_out, *Bb3_out;

	/* Local variables */
	PyObject *div_obj, *phi_obj, *args_obj, *kwargs_obj, \
		 *Bb_obj, *Bb1_arr, *Bb2_arr, *Bb3_arr, *Bb3_tmp, \
	         *Bb1_irrot_arr, *Bb2_irrot_arr, *Bb3_irrot_arr, \
		 *tmp;
	int x0alloc=0;
	npy_intp *dims = NULL, *dimsC = NULL;

	/* Parse input arguments */
	if (!PyArg_ParseTupleAndKeywords(args, kwargs, \
				"OOO|(ddd)(iii)OOdi", killdiv_kwlist, \
				&Bb1_obj, &Bb2_obj, &Bb3_obj, \
				&hx, &hy, &hz, &bX, &bY, &bZ, &boundary, \
				&x0_obj, &thr, &maxiter)) {
		return NULL;
	}

	/* Interprete input objects as `numpy` arrays */
	Bb1_arr = PyArray_FROM_OTF(Bb1_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb2_arr = PyArray_FROM_OTF(Bb2_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);
	Bb3_arr = PyArray_FROM_OTF(Bb3_obj, NPY_DOUBLE, NPY_ARRAY_IN_FARRAY);

	/* Throw an exception if that didn't work */
	if (Bb1_arr == NULL || Bb2_arr == NULL || Bb3_arr == NULL) {
		Py_XDECREF(Bb1_arr);
		Py_XDECREF(Bb2_arr);
		Py_XDECREF(Bb3_arr);
		PyErr_SetString(PyExc_ValueError, \
				"divcleaner_killdiv: Failed parsing `Bb1` \
				and/or `Bb2` and/or `Bb3`.");
		return NULL;
	}

	/* Compute div(B) */
	args_obj    = Py_BuildValue("OOO", Bb1_obj, Bb2_obj, Bb3_obj);
	kwargs_obj  = Py_BuildValue("{s(ddd)}", \
			"grid_spacing", hx, hy, hz);
	if (kwargs_obj == NULL) printf("null...\n");
	div_obj     = divcleaner_div(NULL, args_obj, kwargs_obj);
	Py_XDECREF(args_obj);
	Py_XDECREF(kwargs_obj);
	if (div_obj == NULL) return NULL;

	/* Store array dimensions */
	dims = PyArray_DIMS((PyArrayObject*)Bb1_arr);
	dimsC = PyMem_Malloc(3*sizeof(npy_intp));
	dimsC[0] = dims[0]-1;
	dimsC[1] = dims[1];
	dimsC[2] = dims[2];

	/* Build x0 if it was not provided */
	if (x0_obj == Py_None) {
		tmp = PyArray_ZEROS(3,dimsC,NPY_DOUBLE,1);
		x0_obj = PyNumber_Add(PyFloat_FromDouble(1.),tmp);
		Py_DECREF(tmp);
		x0alloc = 1;
	}

	/* Solve the poisson equation */
	args_obj = Py_BuildValue("O(iii)di", x0_obj, bX, bY, bZ, thr, maxiter);
	kwargs_obj = Py_BuildValue("{sOsOs(ddd)}", "b", div_obj, \
			"boundary", boundary, "grid_spacing", hx, hy, hz);
	phi_obj = divcleaner_poisson(NULL, args_obj, kwargs_obj);
	Py_XDECREF(args_obj);
	Py_XDECREF(kwargs_obj);
	Py_XDECREF(div_obj);
	if (x0alloc) Py_XDECREF(x0_obj);
	if (phi_obj == NULL) return NULL;
	
	/* Compute irrotational field component */
	args_obj = Py_BuildValue("(O)", phi_obj);
	kwargs_obj  = Py_BuildValue("{s(ddd)}", \
			"grid_spacing", hx, hy, hz);
	Bb_obj = divcleaner_phiToBbc(NULL, args_obj, kwargs_obj);
	Py_XDECREF(args_obj);
	Py_XDECREF(kwargs_obj);
	Py_XDECREF(phi_obj);
	if (Bb_obj == NULL) return NULL;
	Bb1_irrot_arr = PyTuple_GetItem(Bb_obj, 0);
	Bb2_irrot_arr = PyTuple_GetItem(Bb_obj, 1);
	Bb3_irrot_arr = PyTuple_GetItem(Bb_obj, 2);
	Py_XDECREF(Bb_obj);

	Bb1_out = (PyObject*)PyNumber_Subtract(Bb1_arr, Bb1_irrot_arr);
	Bb2_out = (PyObject*)PyNumber_Subtract(Bb2_arr, Bb2_irrot_arr);
	Bb3_tmp = (PyObject*)PyNumber_Subtract(Bb3_arr, Bb3_irrot_arr);
	Py_XDECREF(Bb1_arr);
	Py_XDECREF(Bb2_arr);
	Py_XDECREF(Bb3_arr);
	Py_XDECREF(Bb1_irrot_arr);
	Py_XDECREF(Bb2_irrot_arr);
	Py_XDECREF(Bb3_irrot_arr);

	/* Fix top and bottom boundaries */
	args_obj = Py_BuildValue("OOO", Bb1_out, Bb2_out, Bb3_tmp);
	kwargs_obj  = Py_BuildValue("{s(ddd)}", \
			"grid_spacing", hx, hy, hz);
	Bb3_out = divcleaner_fixTopBottomBoundary(NULL, args_obj, kwargs_obj);
	Py_XDECREF(Bb3_tmp);
	Py_XDECREF(args_obj);
	Py_XDECREF(kwargs_obj);
	if (Bb3_out == NULL) return NULL;

	return Py_BuildValue("OOO", Bb1_out, Bb2_out, Bb3_out);
}
