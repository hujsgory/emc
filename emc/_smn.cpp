#ifdef __cplusplus
extern "C" {
#endif

    #include <Python.h>
    #include "numpy/arrayobject.h"

#ifdef __cplusplus
}  // extern "C"
#endif

#include <math.h>
#define PI (3.14159265358979323846)


static inline npy_double sint(npy_double dx,npy_double dy){
    return dy/hypot(dx,dy);
}
static inline npy_double cost(npy_double dx,npy_double dy){
    return dx/hypot(dx,dy);
}
static inline npy_double center(npy_double beg, npy_double end){
    return (beg+end)/2.0;
}
static inline npy_double getBegSubint(npy_double beg,npy_double end,int i,int n){
    return beg+npy_double(i)  *(end-beg)/npy_double(n);
}
static inline npy_double getEndSubint(npy_double beg,npy_double end,int i,int n){
    return beg+npy_double(i+1)*(end-beg)/npy_double(n);
}
static npy_double sumatan(npy_double a1, npy_double a2, npy_double c){
    if (c==0.0)
        return -(a1+a2)/(a1*a2);
    npy_double c2=c*c;
    npy_double a12=a1*a2;
    npy_double atg=0;
    if (c2==a12){
        if (c*a1>0.0)
            atg=PI/2;
        else 
            atg=-PI/2;
    }
    else{
        atg=atan(c*(a1+a2)/(c2-a1*a2));
        if (a12/c2 > 1.0){
            if (c*a1 > 0.0) 
                atg += PI;
            else 
                atg-= PI;
        }
    }
    return atg;
}
static inline npy_double F1(npy_double a, npy_double b, npy_double dn2){
   npy_double k=dn2+b;
   npy_double l=dn2-b;
   if (a!=0.0){
       npy_double k1=k*k+a*a;
       npy_double l1=l*l+a*a;
       return dn2*(log(k1*l1)-4.)+ 2.*a*sumatan(k,l,a)+ b*log(k1/l1);
   }
   return b*log(k*k/(l*l))+dn2*(log(l*l*k*k)-4.);
}
static inline npy_double F2(npy_double a, npy_double b, npy_double dn2){
    if (a==0.0)
        return sumatan(dn2+b, dn2-b, a);
    return     sumatan(dn2+b, dn2-b, a)/a;
}
static inline npy_double F3(npy_double a, npy_double b, npy_double dn2){
    npy_double k=dn2+b;
    npy_double l=dn2-b;
    npy_double ka=k*k+a*a;
    npy_double la=l*l+a*a;
    return 0.5*log(fabs(la/ka));
}
static inline npy_double FI(npy_double a1, npy_double a2, npy_double c, npy_double dn){
    if (c!=0.0)
        return 2*(c*(atan(a1/c)-atan(a2/c))-dn)+a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c);
    return a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)-2*dn;
}

static char _smn_any_doc[]="Value of matrix S calculation";

//self,block_S,list1,list2,bDiel
static PyObject * _smn_any(PyObject * self, PyObject * args){
    npy_bool bDiel=0,iflg=1;
    PyObject  *py_list1=NULL, *py_list2=NULL;
    PyArrayObject *block_S=NULL;
    PyArg_ParseTuple(args,"OOObb", &block_S, &py_list1, &py_list2, &bDiel, &iflg);
    npy_double * data=(npy_double *)PyArray_DATA(block_S);

    Py_ssize_t len_list1 = PyList_Size(py_list1);
    Py_ssize_t len_list2 = PyList_Size(py_list2);
    for (Py_ssize_t idx_list1=0; idx_list1<len_list1; idx_list1++){
        PyObject * bound_m=PyList_GetItem(py_list1,idx_list1);
        PyObject * section_m=PyDict_GetItemString(bound_m,"_section_");
        npy_double m_begx = PyFloat_AsDouble(PyList_GetItem(section_m,0));
        npy_double m_begy = PyFloat_AsDouble(PyList_GetItem(section_m,1));
        npy_double m_endx = PyFloat_AsDouble(PyList_GetItem(section_m,2));
        npy_double m_endy = PyFloat_AsDouble(PyList_GetItem(section_m,3));
        short len_bound_m = (short)PyLong_AsLong((PyObject *)PyDict_GetItemString(bound_m,"n_subint"));
        npy_double sinm   = sint(m_endx-m_begx,m_endy-m_begy);
        npy_double cosm   = cost(m_endx-m_begx,m_endy-m_begy);
        for(short i=0; i<len_bound_m; i++){
            npy_double subs_i_begx=getBegSubint(m_begx,m_endx,i,len_bound_m);
            npy_double subs_i_endx=getEndSubint(m_begx,m_endx,i,len_bound_m);
            npy_double subs_i_begy=getBegSubint(m_begy,m_endy,i,len_bound_m);
            npy_double subs_i_endy=getEndSubint(m_begy,m_endy,i,len_bound_m);
            npy_double xm=center(subs_i_begx,subs_i_endx),ym=center(subs_i_begy,subs_i_endy);
            int n=0;
            for (Py_ssize_t idx_list2=0; idx_list2<len_list2; idx_list2++){
                PyObject * bound_n=PyList_GetItem(py_list2,idx_list2);
                PyObject * section_n=PyDict_GetItemString(bound_n,"_section_");
                npy_double n_begx = PyFloat_AsDouble(PyList_GetItem(section_n,0));
                npy_double n_begy = PyFloat_AsDouble(PyList_GetItem(section_n,1));
                npy_double n_endx = PyFloat_AsDouble(PyList_GetItem(section_n,2));
                npy_double n_endy = PyFloat_AsDouble(PyList_GetItem(section_n,3));
                short len_bound_n = (short)PyLong_AsLong(PyDict_GetItemString(bound_n,"n_subint"));
                npy_double sinn   = sint(n_endx-n_begx,n_endy-n_begy);
                npy_double cosn   = cost(n_endx-n_begx,n_endy-n_begy);
                for (short j=0;j<len_bound_n;j++){
                    npy_double subs_j_begx=getBegSubint(n_begx,n_endx,j,len_bound_n);
                    npy_double subs_j_endx=getEndSubint(n_begx,n_endx,j,len_bound_n);
                    npy_double subs_j_begy=getBegSubint(n_begy,n_endy,j,len_bound_n);
                    npy_double subs_j_endy=getEndSubint(n_begy,n_endy,j,len_bound_n);
                    npy_double xn  = center(subs_j_begx,subs_j_endx),yn=center(subs_j_begy,subs_j_endy);
                    npy_double dn2 = hypot(subs_j_endx-subs_j_begx,subs_j_endy-subs_j_begy)/2.0;
                    npy_double dx  = xm - xn;
                    npy_double a2  = dx*sinn-(ym-yn)*cosn;
                    npy_double b2  = dx*cosn+(ym-yn)*sinn;
                    npy_double a1  = dx*sinn+(ym+yn)*cosn;
                    npy_double b1  = dx*cosn-(ym+yn)*sinn;
                    if (!bDiel){
                        *data    = -F1(a2, b2, dn2);
                        if (iflg)
                          *data +=  F1(a1, b1, dn2);
                    }
                    else{
                        npy_double f2= F2(a2, b2, dn2);
                        npy_double f3= F3(a2, b2, dn2);
                        npy_double Imn= ((dx-b2*cosn)*f2-cosn*f3)*sinm - ((ym-yn-b2*sinn)*f2-sinn*f3)*cosm;
                        if (iflg){
                            npy_double f2 = F2(a1, b1, dn2);
                            npy_double f3 = F3(a1, b1, dn2);
                            Imn += -((dx-b1*cosn)*f2-cosn*f3)*sinm + ((ym+yn+b1*sinn)*f2-sinn*f3)*cosm;
                        }
                        *data=Imn;
                    }
                    data++;
                }
            }
        }
    }
    return Py_BuildValue("i",1);
}

static char _smn_ortho_doc[]="Value of matrix S calculation for orthogonal structure";

static PyObject * _smn_ortho(PyObject * self, PyObject * args){
    return Py_BuildValue("i",1);
}

static PyMethodDef _methods[]={{"any", _smn_any, METH_VARARGS, _smn_any_doc},
                               {"ortho", _smn_ortho, METH_VARARGS, _smn_ortho_doc},
                               {NULL   , NULL,            0,     NULL}};
#ifdef __cplusplus
extern "C" {
#endif

__declspec(dllexport) void init_smn(){
    Py_InitModule("_smn",_methods);
}

#ifdef __cplusplus
}  // extern "C"
#endif
