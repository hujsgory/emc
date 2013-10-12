#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#define PI (3.14159265358979323846)

double sumatan(double a1, double a2, double c){
    if (c==0.0)
        return -(a1+a2)/(a1*a2);
    double c2=c*c;
    double a12=a1*a2;
    double atg=0;
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
double F1(double a, double b, double dn2){
   double k=dn2+b;
   double l=dn2-b;
   double k1,l1;
   if (a!=0.0){
       k1=k*k+a*a;
       l1=l*l+a*a;
       return dn2*(log(k1*l1)-4.)+ 2.*a*sumatan(k,l,a)+ b*log(k1/l1);
   }
   return b*log(k*k/(l*l))+dn2*(log(l*l*k*k)-4.);
}
double F2(double a, double b, double dn2){
    if (a==0.0)
        return sumatan(dn2+b, dn2-b, a);
    return     sumatan(dn2+b, dn2-b, a)/a;
}
double F3(double a, double b, double dn2){
    double k=dn2+b;
    double l=dn2-b;
    double ka=k*k+a*a;
    double la=l*l+a*a;
    return 0.5*log(fabs(la/ka));
}
double FI(double a1, double a2, double c, double dn){
    if (c!=0.0)
        return 2*(c*(atan(a1/c)-atan(a2/c))-dn)+a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c);
    return a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)-2*dn;
}

static char _smn_doc[]="Value of matrix S calculation";
//self,block_S,list1,list2,bDiel
static PyObject * _smn(PyObject * self, PyObject * args){
    char bDiel=0;
    PyListObject  *list1=NULL, *list2=NULL;
    PyArrayObject *block_S;
    PyArg_ParseTuple(args,"000b", &block_S, &list1, &list2, &bDiel);
    int m=0;
    Py_ssize_t len_list1 = PyList_GetSize(list1);
    Py_ssize_t len_list2 = PyList_GetSize(list2);
    for (Py_ssize_t idx_list1=0; idx_list1<len_list1; idx_list1++){
        PyDictObject * bound_m=PyList_GetItem(list1,idx_list1);
        //FIXME: Need convert Section to List
        PyListObject * section_m=PyDict_GetItemString(bound_m,"section");
        PyList_GetItem(section_m,0);
        PyList_GetItem(section_m,1);
        PyList_GetItem(section_m,2);
        PyList_GetItem(section_m,3);
    //    sinm,cosm=section_m.sint,section_m.cost
    //    len_bound_m=bound_m['n_subint']
    //    for i in xrange(len_bound_m):
    //        subsection_i=section_m.getSubinterval(i,len_bound_m)
    //        xm,ym = subsection_i.center.x,subsection_i.center.y
            int n=0;
            for (Py_ssize_t idx_list1=0; idx_list1<len_list1; idx_list1++){
    //            section_n = bound_n['section']
    //            sinn,cosn=section_n.sint,section_n.cost
    //            len_bound_n=bound_n['n_subint']
                for (int j=0;j<len_bound_n;j++){
                    //subsection_j=section_n.getSubinterval(j,len_bound_n)
                    double dn2 = subsection_j.len/2;
                    //double xn=subsection_j.center.x,yn=subsection_j.center.y;
                    double dx = xm - xn;
                    double a2=dx*sinn-(ym-yn)*cosn;
                    double b2=dx*cosn+(ym-yn)*sinn;
                    double a1= dx*sinn+(ym+yn)*cosn;
                    double b1= dx*cosn-(ym+yn)*sinn;
                    if (!bDiel){
                        //block_S[m, n] = -self.F1(a2, b2, dn2)
                        //if self.iflg:
                        //  block_S[m, n] += self.F1(a1, b1, dn2)
                    }
                    else{
                        double f2= F2(a2, b2, dn2);
                        double f3= F3(a2, b2, dn2);
                        double Imn= ((dx-b2*cosn)*f2-cosn*f3)*sinm - ((ym-yn-b2*sinn)*f2-sinn*f3)*cosm;
                        if (iflg){
                            double f2 = F2(a1, b1, dn2);
                            double f3 = F3(a1, b1, dn2);
                            Imn += -((dx-b1*cosn)*f2-cosn*f3)*sinm + ((ym+yn+b1*sinn)*f2-sinn*f3)*cosm;
                        }
                        //block_S[m, n]= Imn
                    }
                    n+=1;
                }
            m+=1;
        }
    }

}

static PyMethodDef _smnmethods[]={{"_calcSmn_", _smn, METH_VARARGS, _smn_doc},
                                  {NULL   , NULL,            0,     NULL}};
void init_smn(){
    Py_InitModule("_smn",_smnmethods);
}
