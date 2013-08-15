#!/usr/bin/python
#coding: utf8
from math import *

class Coord(object):
    def __init__(self,x=0.0,y=0.0):
        if type(x) and type(y) is float:
            self.x,self.y=x,y
        else: raise TypeError
    def __eq__(self,_coord):
        return self.x==_coord.x and self.y==_coord.y
    def __ne__(self,_coord):
        return self.x!=_coord.x or self.y!=_coord.y
    def sint(self,_coord):
        return (_coord.y-self.y)/hypot(_coord.y-self.y,_coord.x-self.x)
    def cost(self,_coord):
        return (_coord.x-self.x)/hypot(_coord.y-self.y,_coord.x-self.x)
    def center(self,_coord):
        return Coord((_coord.x+self.x)/2.0,(_coord.y+self.y)/2.0)
    def len(self,_coord):
        return hypot(_coord.y-self.y,_coord.x-self.x)

class Section(object):
    def __init__(self,beg=Coord(),end=Coord()):
        if type(beg) and type(end) is Coord:
            self.beg,self.end=beg,end
        else: raise TypeError
    @property
    def center(self):
        return self.beg.center(self.end)
    @property
    def sint(self):
        return self.beg.sint(self.end)
    @property
    def cost(self):
        return self.beg.cost(self.end)
    @property
    def len(self):
        return self.beg.len(self.end)
    @property
    def dx(self):
        return self.end.x-self.beg.x
    @property
    def dy(self):
        return self.end.y-self.beg.y
    
class Conf(object):
    def __init__(self):
        self.ListBound=[]
        self.iflg=True
    def __iter__(self):
        return self.ListBound.__iter__()
    def intersection(self,sect2):
        for bound in self.ListBound:
            sect1=bound['section']
            # Coeficients of Ax+By+D=0
            a1=-sect1.dy;
            b1= sect1.dx;
            d1=-(a1*sect1.beg.x+b1*sect1.beg.y);
            a2=-sect2.dy;
            b2= sect2.dx;
            d2=-(a2*sect2.beg.x+b2*sect2.beg.y);

            # Calculate halfplane of sections ends
            seg1_line2_start = a2*sect1.beg.x + b2*sect1.beg.y + d2;
            seg1_line2_end   = a2*sect1.end.x + b2*sect1.end.y + d2;
            seg2_line1_start = a1*sect2.beg.x + b1*sect2.beg.y + d1;
            seg2_line1_end   = a1*sect2.end.x + b1*sect2.end.y + d1;
            
            h1=seg1_line2_start*seg1_line2_end;
            h2=seg2_line1_start*seg2_line1_end;

            # Ends located in the different halfplane
            if h1<0. and h2<0.: return True
            # Collinear (both ends of segments lie on one line)
            if a1==a2 and b1==b2 and d1==d2:
                #   |------|======|------|
                # fmin1  fmin2  fmax1  fmax2
                fmin1=min(sect1.beg.x,sect1.end.x)
                fmin2=min(sect2.beg.x,sect2.end.x)
                fmax1=max(sect1.beg.x,sect1.end.x)
                fmax2=max(sect2.beg.x,sect2.end.x)
                if fmin1<fmax2 and fmin2<fmax1:
                    return True
                fmin1=min(sect1.beg.y,sect1.end.y)
                fmin2=min(sect2.beg.y,sect2.end.y)
                fmax1=max(sect1.beg.y,sect1.end.y)
                fmax2=max(sect2.beg.y,sect2.end.y)
                if fmin1<fmax2 and fmin2<fmax1:
                    return True
        return False
    #
    # @param mat_type: False - Conductor-Dielectric bound, True - Dielectric-Dielectric bound
    # @param n_subint: number of bound's subintervals
    # @param mat_param: erp - relative permittivity on right side of section, erm - on left side; tdp,tdm - tangent dielectric loss; 
    def add(self,section,mat_type,n_subint=1,**mat_param):
        if (type(section) is Section) and (type(mat_type) is bool) and (type(n_subint) is int) and (type(mat_param) is dict):
            if self.intersection(section): raise ValueError
            erp = mat_param.get('erp', 1.0)
            erm = mat_param.get('erm', 1.0)
            if erp!=erm:
                self.ListBound.append({'section':section,'mat_type':mat_type,'n_subint':n_subint,'mat_param':mat_param})
            else: raise ValueError
        else: raise TypeError



def a1(m,n):
    if type(m) and type(n) is Section:
        return (m.center.x-n.center.x)*n.sint+(m.center.y+n.center.y)*n.cost
    else: raise TypeError
def a2(m,n):
    if type(m) and type(n) is Section:
        return (m.center.x-n.center.x)*n.sint-(m.center.y-n.center.y)*n.cost
    else: raise TypeError
def b1(m,n):
    if type(m) and type(n) is Section:
        return (m.center.x-n.center.x)*n.cost-(m.center.y+n.center.y)*n.sint
    else: raise TypeError
def b2(m,n):
    if type(m) and type(n) is Section:
        return (m.center.x-n.center.x)*n.cost+(m.center.y-n.center.y)*n.sint
    else: raise TypeError

'''
def F1(a,b,dn):
    if type(a) and type(b) and type(dn) is float:
        t1=dn+2.0*b
        t2=dn-2.0*b
        t3=2.0*a
        t4=t1**2.0+t3**2.0
        t5=t2**2.0+t3**2.0
        return b*log(t4/t5)+dn/2.0*log(t4*t5/16.0)-2.0*dn+t3*(atan2(t2,t3)+atan2(t1,t3))
    else:raise TypeError

def F2(a,b,dn):
    if type(a) and type(b) and type(dn) is float:
        if a==0.0:
            return 4.0*dn/(4.0*b**2.0+dn**2.0)
        else:
            return (atan2(dn-2.0*b,2.0*a)+atan2(dn+2.0*b,2.0*a))/a
    else:raise TypeError

def F3(a,b,dn):
    if type(a) and type(b) and type(dn) is float:
        return log(((dn-2.0*b)**2.0+4.0*a**2.0)/((dn+2.0*b)**2.0+4.0*a**2.0))/2.0
    else:raise TypeError
'''
def Imn(m,n):
    if type(m) and type(n) is Section:
        return (m.center.y-n.center.y-b2(m,n)*n.sint)*F2(a2(m,n),b2(m,n),n.len)-n.sint*F3(a2(m,n),b2(m,n),n.len)
    else:raise TypeError

def I_mn(m,n):
    if type(m) and type(n) is Section:
        return (m.center.y+n.center.y+b1(m,n)*n.sint)*F2(a1(m,n),b1(m,n),n.len)-n.sint*F3(a1(m,n),b1(m,n),n.len)
    else:raise TypeError
def Smn(config):
    if type(config) is Conf:
        for bound in config:
            erp = bound['mat_param'].get('erp', 1.0)
            erm = bound['mat_param'].get('erm', 1.0)
            tdp = bound['mat_param'].get('tdp', 0.0)
            tdm = bound['mat_param'].get('tdm', 0.0)
        #return []
    else: raise TypeError

def sumatan(a1,a2,c):
	if c==0.0:
		return -(a1+a2)/(a1*a2)
	c2=c*c
	a12=a1*a2
	if c2==a12:
		if c*a1>0.0:
			atg=pi/2
		else:
			atg=-pi/2
	else:
		atg=atan(c*(a1+a2)/(c2-a1*a2))
		if a12/c2 > 1.:
			if c*a1 > 0.0:
				atg += pi
			else:
				atg-= pi
	return atg

def F1(a,b,dn2):
	k=dn2+b
	l=dn2-b
	if a!=0.0:
		k1=k*k+a*a
		l1=l*l+a*a
		return dn2*(log(k1*l1)-4.)+ 2.*a*sumatan(k,l,a)+ b*log(k1/l1)
	return b*log(k*k/(l*l))+dn2*(log(l*l*k*k)-4.)
def F2(a,b,dn2):
	if a==0.0:
		return sumatan(dn2+b, dn2-b, a)
	return  sumatan(dn2+b,dn2-b,a)/a

def F3(a,b,dn2):
	k=dn2+b
	l=dn2-b
	k= k*k+a*a
	l= l*l+a*a
	return 0.5*log(fabs(l/k))

def FI(a1, a2, c, dn):
	if c!=0.0:
		return 2*(c*(atan(a1/c)-atan(a2/c))-dn)+a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)
	return a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)-2*dn
