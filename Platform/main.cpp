#include "Python.h"
#include <cstdio>

int main(int argc, char **argv)
{
    Py_Initialize();
    PySys_SetArgv(argc, argv);
    PyObject* PyFileObject = PyFile_FromString("main.py", "r");
    if (PyFileObject>0){
        PyRun_SimpleFileEx(PyFile_AsFile(PyFileObject), "0", 1);
    }
    Py_Finalize();
    return 0;
}

