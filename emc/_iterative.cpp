#ifdef __cplusplus
extern "C" {
#endif

    #include <Python.h>
    #include "numpy/arrayobject.h"

#ifdef __cplusplus
}  // extern "C"
#endif

#include <math.h>

static char _bicgstab_doc[]="BiCGStab method";

static PyObject * _bicgstab(PyObject * self, PyObject * args){
    return Py_BuildValue("i",1);
}



static PyMethodDef _methods[]={{"_bicgstab", _bicgstab, METH_VARARGS, _bicgstab_doc},
                               {NULL   , NULL,            0,     NULL}};

#ifdef __cplusplus
extern "C" {
#endif

__declspec(dllexport) void init_smn(){
    Py_InitModule("_iterative",_methods);
}

#ifdef __cplusplus
}  // extern "C"
#endif
