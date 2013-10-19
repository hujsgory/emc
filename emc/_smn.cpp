#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#define PI (3.14159265358979323846)

static inline double hypot(double dx, double dy){
	return sqrt(dx*dx+dy*dy);
}
static inline double sint(double dx,double dy){
	return dx/hypot(dx,dy);
}
static inline double sint(double dx,double dy){
	return dy/hypot(dx,dy);
}
static inline double center(double beg, double end){
	return (beg+end)/2.0;
}
static inline double getBegSubint(double beg,double end,int i,int n){
	return beg+double(i)  *(end-beg)/double(n);
}
static inline double getEndSubint(double beg,double end,int i,int n){
	return beg+double(i+1)*(end-beg)/double(n);
}
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
        PyListObject * section_m=PyDict_GetItemString(bound_m,"_section_");
        double m_begx=PyFloat_AsDouble(PyList_GetItem(section_m,0));
        double m_begy=PyFloat_AsDouble(PyList_GetItem(section_m,1));
        double m_endx=PyFloat_AsDouble(PyList_GetItem(section_m,2));
        double m_endy=PyFloat_AsDouble(PyList_GetItem(section_m,3));
        double sinm = sint(m_endx-m_begx,m_endy-m_begy);
        double cosm = cost(m_endx-m_begx,m_endy-m_begy);
        int len_bound_m=(int)PyLong_AsLong(PyDict_GetItemString(bound_m,"n_subint");
        for(int i=0; i<len_bound_m; i++){
            double subs_i_begx=getBegSubint(m_begx,m_endx,i,len_bound_m);
			double subs_i_endx=getEndSubint(m_begx,m_endx,i,len_bound_m);
			double subs_i_begy=getBegSubint(m_begy,m_endy,i,len_bound_m);
			double subs_i_endy=getEndSubint(m_begy,m_endy,i,len_bound_m);
            double xm=center(subs_i_begx,subs_i_endx),ym=center(subs_i_begy,subs_i_endy);
            int n=0;
            for (Py_ssize_t idx_list2=0; idx_list2<len_list2; idx_list2++){
				PyDictObject * bound_n=PyList_GetItem(list2,idx_list2);
				PyListObject * section_n=PyDict_GetItemString(bound_n,"_section_");
				double n_begx=PyFloat_AsDouble(PyList_GetItem(section_n,0));
				double n_begy=PyFloat_AsDouble(PyList_GetItem(section_n,1));
				double n_endx=PyFloat_AsDouble(PyList_GetItem(section_n,2));
				double n_endy=PyFloat_AsDouble(PyList_GetItem(section_n,3));
				double sinn = sint(n_endx-n_begx,n_endy-n_begy);
				double cosn = cost(n_endx-n_begx,n_endy-n_begy);
				int len_bound_n=(int)PyLong_AsLong(PyDict_GetItemString(bound_n,"n_subint");
                for (int j=0;j<len_bound_n;j++){
					double subs_j_begx=getBegSubint(n_begx,n_endx,j,len_bound_n);
					double subs_j_endx=getEndSubint(n_begx,n_endx,j,len_bound_n);
					double subs_j_begy=getBegSubint(n_begy,n_endy,j,len_bound_n);
					double subs_j_endy=getEndSubint(n_begy,n_endy,j,len_bound_n);
					double xm=center(subs_j_begx,subs_j_endx),ym=center(subs_j_begy,subs_j_endy);
                    double dn2 = hypot(subs_j_endx-subs_j_begx,subs_j_endy-subs_j_begy)/2.0;
                    double dx = xm - xn;
                    double a2=dx*sinn-(ym-yn)*cosn;
                    double b2=dx*cosn+(ym-yn)*sinn;
                    double a1= dx*sinn+(ym+yn)*cosn;
                    double b1= dx*cosn-(ym+yn)*sinn;
                    if (!bDiel){
                        //block_S[m, n] = -F1(a2, b2, dn2);
                        //if (iflg)
                        //  block_S[m, n] += F1(a1, b1, dn2);
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
