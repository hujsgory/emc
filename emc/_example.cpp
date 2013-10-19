#ifdef __cplusplus
extern "C" {
#endif 
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdio.h>
#ifdef __cplusplus
}  // extern "C"
#endif
static char _calc_doc[]="Value of matrix entries calculation";

static PyObject * _calc(PyObject * self, PyObject * args){
	PyArrayObject *block_S;
	PyArg_ParseTuple(args,"O", &block_S);
	int ndim=PyArray_NDIM(block_S);
	npy_intp * dims=PyArray_DIMS(block_S);
	int len_array=1;
	for (int i=0; i<ndim; i++){
		len_array*=dims[i];
	}
	double * data=(double *)PyArray_DATA(block_S);
	
	for (int i=0; i<len_array;i++){
		data[i]=0.0;
	}
	return Py_BuildValue("i",ndim);
}

static PyMethodDef _examplemethods[]={{"_calc", _calc, METH_VARARGS, _calc_doc},
									  {NULL    , NULL,             0,     NULL}};

#ifdef __cplusplus
extern "C" {
#endif 
__declspec(dllexport) void init_example(){
    Py_InitModule("_example",_examplemethods);
}


#ifdef __cplusplus
}  // extern "C"
#endif
