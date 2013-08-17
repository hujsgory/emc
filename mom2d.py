#!/usr/bin/python
#coding: utf8
from math import *
import numpy

class Coord(object):
    def __init__(self,x=0.0,y=0.0):
        if type(x) is float and type(y) is float:
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
        if type(beg) is Coord and type(end) is Coord:
            self.beg,self.end=beg,end
        else: raise TypeError
    def __eq__(self,_section):
        return self.beg==_section.beg and self.end==_section.end
    def __ne__(self,_section):
        return self.beg!=_section.beg or self.end!=_section.end
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
    def getSubinterval(self,i=0,n=1):
        if type(i) is int and type(n) is int:
            if i>=n or i<0: raise ValueError
            x,y=self.beg.x,self.beg.y
            dx,dy=self.dx,self.dy
            return Section(Coord(x+i*dx/n,y+i*dy/n),Coord(x+(i+1)*dx/n,y+(i+1)*dy/n))
        else: raise TypeError

class Conf(object):
    def __init__(self):
        self.listBounds=[]
        self.iflg=True
    def __iter__(self):
        return self.listBounds.__iter__()
    def intersection(self,sect2):
        for bound_m in self.listBounds:
            sect1=bound_m['section']
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
    # @param mat_type: False - Conductor-Dielectric bound_m, True - Dielectric-Dielectric bound_m
    # @param n_subint: number of bound_m's subintervals
    # @param mat_param: erp - relative permittivity on right side of section, erm - on left side; tdp,tdm - tangent dielectric loss; 
    def add(self,section,mat_type,n_subint=1,**mat_param):
        if type(section) is Section and type(mat_type) is bool and type(n_subint) is int and type(mat_param) is dict:
            if self.intersection(section): raise ValueError
            erp = mat_param.get('erp', 1.0)
            erm = mat_param.get('erm', 1.0)
            if erp!=erm:
                self.listBounds.append({'section':section,'mat_type':mat_type,'n_subint':n_subint,'mat_param':mat_param})
            else: raise ValueError
        else: raise TypeError

'''
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
        for bound_m in config:
            erp = bound_m['mat_param'].get('erp', 1.0)
            erm = bound_m['mat_param'].get('erm', 1.0)
            tdp = bound_m['mat_param'].get('tdp', 0.0)
            tdm = bound_m['mat_param'].get('tdm', 0.0)
        #return []
    else: raise TypeError
'''

'''
Port from smn.cpp
'''

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

def SmnAny2D(conf):
    if type(conf) is not Conf:
        raise TypeError
    bUpdate=False
    conf.listBounds.sort(key=lambda x: x['mat_type'])
    m_size = reduce(lambda r,x: r+x['n_subint'],conf.listBounds,0)
    if not conf.iflg :
        m_size+=1 # add one row and one column
    ret=numpy.zeros((m_size,m_size))
    m = 0
    for bound_m in conf.listBounds:
        section_m=bound_m['section']
        n_subint_m=bound_m['n_subint']
        # check if conductors was processed
        bDiel = bound_m['mat_type']
        # set parameters equal for all subsections
        sinm,cosm=section_m.sint,section_m.cost
        dxm,dym=section_m.dx/n_subint_m, section_m.dy/n_subint_m
        first_subsection=section_m.getSubinterval(0,n_subint_m)
        # xm, ym - centers of subinterval
        xm,ym = first_subsection.center.x, first_subsection.center.y
        er_plus=0.0
        erp = bound_m['mat_param'].get('erp', 1.0)
        erm = bound_m['mat_param'].get('erm', 1.0)
        if bDiel : # edit Er for dielectrics
            # Pascal: p -> u, m -> d
            if erp==erm:
                raise ValueError
            er_plus=(erp+erm)*pi/(erp-erm)
        # BEGIN cycle through subintervals
        for subint_m in xrange(n_subint_m):
            # DO JUST THE SAME CALCULATIONS FOR INTEGRAL INTERVALS
            n = 0
            for bound_n in conf.listBounds:
                section_n = bound_n['section'] # BEGIN cycle through integral intervals
                n_subint_n=bound_n['n_subint']
                # set parameters equal for all subsections
                sinn,cosn= section_n.sint, section_n.cost
                dxn,dyn=section_n.dx/n_subint_n, section_n.dy/n_subint_n
                dn = section_n.len/n_subint_n
                dn2 = dn/2.0
                a1,b1=0.0,0.0
                first_subsection_i=section_n.getSubinterval(0,n_subint_n)
                # xn, yn - centers of subinterval
                xn,yn = first_subsection_i.center.x, first_subsection_i.center.y
                dx = xm - xn
                a2 = dx*sinn-(ym-yn)*cosn
                b2 = dx*cosn+(ym-yn)*sinn
                if conf.iflg:
                    a1= dx*sinn+(ym+yn)*cosn
                    b1= dx*cosn-(ym+yn)*sinn
                # BEGIN cycle through integral subintervals
                for subint_n in xrange(n_subint_n):
                    # if updating existing smn matrix, we need to update the diagonal only, skip calculations for non-diagonal elements
                    if not bUpdate or m==n: 
                        if not bDiel:
                            # CALCULATION WITH CONDUCTORS
                            ret[m, n] = -F1(a2, b2, dn2)
                            if conf.iflg:
                                ret[m, n] += F1(a1, b1, dn2)
                                b1 -= dn
                        else :
                            # CALCULATION WITH DIELECTRICS
                            # Imn= sinm*(Ix-Ix1)-cosm*(Iy-Iy1)
                            f2= F2(a2, b2, dn2)
                            f3= F3(a2, b2, dn2)
                            Imn= 0.0
                            if sinm!=0.0:
                                Imn+=((dx   -b2*cosn)*f2-cosn*f3)*sinm  #  Ix*sinm
                            if cosm!=0.0:
                                Imn-=((ym-yn-b2*sinn)*f2-sinn*f3)*cosm  # -Iy*cosm
                            if conf.iflg:
                                f2= F2(a1, b1, dn2)
                                f3= F3(a1, b1, dn2)
                                if sinm!=0.0:
                                    Imn -= ((dx   -b1*cosn)*f2-cosn*f3)*sinm  # -Ix1*sinm
                                if cosm!=0.0:
                                    Imn += ((ym+yn+b1*sinn)*f2-sinn*f3)*cosm  #+Iy1*cosm
                                # Iy1 - I1, Iy - I2, 
                                b1-= dn
                            if m==n:
                                Imn += er_plus
                            if Imn == 0.0:
                                ret[m, n] = -1.0
                            ret[m, n]= Imn
                    else:
                        if conf.iflg :
                            b1-= dn
                    b2-= dn
                    dx-= dxn
                    n+=1 # increment matrix index
                    xn += dxn # calc center of next integral subsection
                    yn += dyn
                 # END cycle through integral subintervals
             # END cycle through integral intervals
            m+=1 # increment matrix index
            xm += dxm # calc center of next subsection
            ym += dym
         # END cycle through subintervals
     # END cycle through intervals
    #   if(!conf->_bInfiniteGround) : # fill in additional row and column
    if not conf.iflg : # fill in additional row and column
        sz = m_size-1
        n = 0
        for bound_m in conf.listBounds: # cycle through conductors
            section_m=bound_m['section']
            n_subint_m=bound_m['n_subint']
            for si in xrange(n_subint_m):
                if bound_m['mat_type']==False:
                    si_len = section_m.len/n_subint_m
                    ret[n, sz]=si_len/ret[n, n]
                    ret[sz, n]=si_len*bound_m['mat_param'].get('erp', 1.0)
                else: # clear rest of matrix cells
                    ret[n,sz],ret[sz,n] = 0.0,0.0
                n+=1
    return ret

def SmnOrtho(conf) :
    if type(conf) is not Conf:
        raise TypeError
    bUpdate=False
    conf.listBounds.sort(key=lambda x: x['mat_type'])
    m_size = reduce(lambda r,x: r+x['n_subint'],conf.listBounds,0)
    if not conf.iflg: 
        ++m_size; # add one row and one column
    ret=numpy.zeros((m_size,m_size))
    m = 0
    for bound_m in conf: # section_mIN cycle through intervals
        section_m=bound_m['section']
        n_subint_m=bound_m['n_subint']
        # check if conductors was processed
        bDiel = bound_m['mat_type']
        # set parameters equal for all subsections
        sinm = section_m.sint
        negcosm = -section_m.cost 
        dxm,dym=section_m.dx/n_subint_m, section_m.dy/n_subint_m
        first_subsection=section_m.getSubinterval(0,n_subint_m)
        er_plus = 0
        # xm, ym - centers of subinterval
        xm,ym = first_subsection.center.x, first_subsection.center.y
        er_plus=0.0
        erp = bound_m['mat_param'].get('erp', 1.0)
        erm = bound_m['mat_param'].get('erm', 1.0)
        if bDiel: # edit Er for dielectrics
            # Pascal: p -> u, m -> d
            if erp==erm:
                raise ValueError
            er_plus=(erp+erm)*pi/(erp-erm)

        # BEGIN cycle through subintervals
        for subint_m in xrange(n_subint_m):
            # DO JUST THE SAME CALCULATIONS FOR INTEGRAL INTERVALS
            n = 0
            for bound_n in conf.listBounds:
                section_n = bound_n['section'] # BEGIN cycle through integral intervals
                n_subint_n=bound_n['n_subint']
                # set parameters equal for all subsections
                sinn = section_n.sint
                cosn = section_n.cost
                dxn,dyn=section_n.dx/n_subint_n, section_n.dy/n_subint_n
                dn = section_n.len/n_subint_n
                dn2 = dn/2.0
                a1 = 0. #/*b1 =0, */
                c1,c2,da,fi=0.,0.,0.,0.
                first_subsection_i=section_n.getSubinterval(0,n_subint_n)
                # xn, yn - centers of subinterval
                xn,yn = first_subsection_i.center.x, first_subsection_i.center.y
                dx = xm - xn
                a2 = dx*sinn - (ym-yn)*cosn
                #/*,b2 = dx*cosn + (ym-yn)*sinn*/
                if sinn == 0.0 : # ortho to Y: I == I1
                    a1= dn2-xm+xn
                    a2= a1-dn
                    c1= ym- yn
                    c2= ym+yn
                    # WONTFIX: Refactoring
                    da= dn * section_n.cost  # +-dn
                else :  # ortho X: I - I1 later
                    a1= dn2-ym+yn  # just for I, later more
                    a2= a1- dn
                    c1= c2= xm-xn
                    da= dn*sinn    # +-dn

                # BEGIN cycle through integral subintervals
                for subint_n in xrange(n_subint_n):
                    if not bUpdate or m==n:
                        if not bDiel:
                        # CALCULATION WITH CONDUCTORS
                            fi= -FI(a1, a2, c1, dn)
                            if conf.iflg :
                                if sinn==0.0 :
                                    fi+=FI(a1, a2, c2, dn)
                                else :
                                    fi+=FI(a1+2*ym, a2+2*ym, c2, dn)
                        else :
                            # CALCULATION WITH DIELECTRICS
                            if sinm==0.0: #  m- ortho Y
                                if sinn==0.0:
                                    # both ortho Y
                                    if c1:
                                        fi= sumatan(a1,-a2,c1)
                                    else:
                                        fi= 0.0
                                    if conf.iflg and c2!=0.0:
                                        fi-= sumatan(a1,-a2,c2)
                                
                                else :  # m- ortho Y, n- ortho X
                                    tmp3= c1*c1; #HERE с1==с2
                                    fi= 0.5*log((a2*a2+tmp3)/(a1*a1+tmp3))
                                    if conf.iflg :
                                        tmp1= a1+2*ym
                                        tmp2= tmp1- dn
                                        fi-= 0.5*log((tmp2*tmp2+tmp3)/(tmp1*tmp1+tmp3))
                                # HACK: subintervals direction fix by T.R. Gazizov - DISABLED
                                fi*= negcosm
                            else : #  m- ortho X
                                if sinn==0.0 : # m- ortho X, n- ortho Y
                                    fi= 0.5*log((a2*a2+c1*c1)/(a1*a1+c1*c1))
                                    if conf.iflg :
                                        fi-= 0.5*log((a2*a2+c2*c2)/(a1*a1+c2*c2))

                                else :  # both ortho X
                                    if c1 :
                                        fi= sumatan(a1,-a2,c1)
                                    else :
                                        fi= 0.0
                                    if conf.iflg and c2!=0.0:
                                        fi-= sumatan(a1+2*ym,-a2-2*ym,c2)
                                # HACK: subintervals direction fix by T.R. Gazizov - DISABLED
                                fi*= sinm
                            if(m==n):
                                fi+= er_plus
                         # END DIELECTRIC PROCESSING
                        
                        ret[m, n] = fi
                    a1 += da
                    a2 += da
                    n+=1 # increment matrix index
                    xn += dxn # calc center of next integral subsection
                    yn += dyn
                 # END cycle through integral subintervals
             # END cycle through integral intervals
            m+=1 # increment matrix index
            xm += dxm; # calc center of next subsection
            ym += dym
         # END cycle through subintervals
     # END cycle through intervals
    if not conf.iflg : # fill in additional row and column
        sz = m_size-1
        n = 0
        for bound_m in conf.listBounds: # cycle through conductors
            section_m=bound_m['section']
            n_subint_m=bound_m['n_subint']
            for si in xrange(n_subint_m):
                if bound_m['mat_type']==False:
                    si_len = section_m.len/n_subint_m
                    ret[n, sz]=si_len/ret[n, n]
                    ret[sz, n]=si_len*bound_m['mat_param'].get('erp', 1.0)
                else: # clear rest of matrix cells
                    ret[n,sz],ret[sz,n] = 0.0,0.0
                n+=1
    return ret

