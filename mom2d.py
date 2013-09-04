#!/usr/bin/python
#coding: utf8
from math import *
import numpy
import numpy.linalg as la

eps0=8.854187817e-12
Coef_C = 4*pi*eps0
V0 = 299792458.0

class Coord(object):
    def __init__(self,x=0.0,y=0.0):
        self.x,self.y=float(x),float(y)
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
        self.list_bounds=list()
        self.iflg=True
        self.mat_type=False     #mat_type: False - Conductor-Dielectric bound_m, True - Dielectric-Dielectric bound_m
        self.mat_count=0
        self.sect_count=0
        self.mat_param=dict()
    def __iter__(self):
        return self.list_bounds.__iter__()
    def intersection(self,sect2):
        for bound_m in self.list_bounds:
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
            if atan2(a1,b1)==atan2(a2,b2) and d1==d2:
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

    # n_subint: number of bound_m's subintervals
    def add(self,section,n_subint=1):
        if type(section) is Section and type(n_subint) is int:
            if self.intersection(section): raise ValueError
            erp = self.mat_param.get('erp', 1.0)
            erm = self.mat_param.get('erm', 1.0)
            if erp!=erm:
                self.list_bounds.append({'section':section,'mat_type':self.mat_type,'n_subint':n_subint,'mat_param':self.mat_param, 'mat_count': self.mat_count,'sect_count': self.sect_count})
                self.sect_count+=1
            else: raise ValueError
        else: raise TypeError
    # mat_param: erp - relative permittivity on right side of section, erm - on left side; tdp,tdm - tangent dielectric loss; 
    def cond(self,**mat_param):
        self.mat_param=mat_param
        self.mat_type=False
        self.mat_count+=1
        self.sect_count=0
    def diel(self,**mat_param):
        self.mat_param=mat_param
        self.mat_type=True
        self.mat_count+=1
        self.sect_count=0


class Board():
    def __init__(self):
        self.layers=list()
        self.environment=1.0
    def layer(self,height,er,td=0.0):
        if height<=0.0 or er<1.0 or td <0.0:
            raise ValueError
        if len(self.layers)>0:
            if height<=max(map(lambda x: x['thickness']-x['depth'], self.layers[-1]['cond'])):
                raise ValueError('Thickness of conductor of previous layer is greater than height')
            if self.layers[-1]['cover']:
                raise ValueError('Layer can\'t be applied to cover')
        self.layers.append({'height':height,'er':er,'td':td,'cover':False,'cond':list()})
    def conductor(self,space,width,thickness,depth=0.0,to_center=False):
        if len(self.layers)==0:
            raise ValueError('It is necessary to create at least one layer')
        last_layer_cond=self.layers[-1]['cond']
        if to_center:
            space-=width/2.0
            if len(last_layer_cond)>0:
                space-=last_layer_cond[-1]['width']/2.0
        if space<=0. or width<=0. or thickness<=0. or depth<0.:
            raise ValueError('All parameters must be positive')
        if depth>=self.layers[-1]['height']:
            raise ValueError('Depth is greater than layer height')
        last_layer_cond.append({'space':space,'width':width,'thickness':thickness,'depth':depth})
    def cover(self,height,er,td=0.0):
        self.layers.append({'height':height,'er':er,'td':td,'cover':False,'cond':list()})
    def board2conf(self):
        self.conf=Conf()
        y_layer=0.0
        # Calculation of structure's right coordinate x
        self.max_x=max(map(lambda x: reduce(lambda r,y: r+y['space']+y['width'],x['cond'],x['cond'][0]['space']),self.layers))

        
'''
Port from smn.cpp
'''

class Smn(object):
    # list_bounds,iflg,matrix_S,m_size
    def __init__(self,conf):
        if type(conf) is not Conf:
            raise TypeError
        self.list_bounds=sorted(conf.list_bounds,key=lambda x: [x['mat_type'],x['mat_count'],x['sect_count']])
        self.iflg=conf.iflg
        self.n_c,self.n_d=0,0
        for bound in self.list_bounds:
            if bound['mat_type']:
                self.n_d+=bound['n_subint']
            else:
                self.n_c+=bound['n_subint']
        self.m_size=self.n_c+self.n_d
        if not self.iflg :
            self.m_size+=1 # add one row and one column

    def sumatan(self,a1,a2,c):
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
    
    def F1(self,a,b,dn2):
        k=dn2+b
        l=dn2-b
        if a!=0.0:
            k1=k*k+a*a
            l1=l*l+a*a
            return dn2*(log(k1*l1)-4.)+ 2.*a*self.sumatan(k,l,a)+ b*log(k1/l1)
        return b*log(k*k/(l*l))+dn2*(log(l*l*k*k)-4.)
    def F2(self,a,b,dn2):
        if a==0.0:
            return self.sumatan(dn2+b, dn2-b, a)
        return  self.sumatan(dn2+b,dn2-b,a)/a
    
    def F3(self,a,b,dn2):
        k=dn2+b
        l=dn2-b
        k= k*k+a*a
        l= l*l+a*a
        return 0.5*log(fabs(l/k))
    
    def FI(self,a1, a2, c, dn):
        if c!=0.0:
            return 2*(c*(atan(a1/c)-atan(a2/c))-dn)+a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)
        return a1*log(a1*a1+c*c)-a2*log(a2*a2+c*c)-2*dn

    # TODO: Refactoring
    def SmnAny2D(self):
        self.matrix_S=numpy.zeros((self.m_size,self.m_size))
        self.diag_S11_C=numpy.zeros((self.n_d))
        self.diag_S11_L=numpy.zeros((self.n_d))
        bUpdate=False
        m = 0
        for bound_m in self.list_bounds:
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
            er_plus,mu_plus=0.0,0.0
            erp = bound_m['mat_param'].get('erp', 1.0)
            erm = bound_m['mat_param'].get('erm', 1.0)
            mup = bound_m['mat_param'].get('mup', 1.0)
            mum = bound_m['mat_param'].get('mum', 1.00001)
            if erp==erm:
                raise ValueError('Dielectric constant of right side is equal to value of left side')
            er_plus=(erp+erm)*pi/(erp-erm)
            mu_plus=(1/mup+1/mum)*pi/(1/mup-1/mum)
            # BEGIN cycle through subintervals
            for subint_m in xrange(n_subint_m):
                # DO JUST THE SAME CALCULATIONS FOR INTEGRAL INTERVALS
                n = 0
                for bound_n in self.list_bounds:
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
                    if self.iflg:
                        a1= dx*sinn+(ym+yn)*cosn
                        b1= dx*cosn-(ym+yn)*sinn
                    # BEGIN cycle through integral subintervals
                    for subint_n in xrange(n_subint_n):
                        if not bDiel:
                            # CALCULATION WITH CONDUCTORS
                            self.matrix_S[m, n] = -self.F1(a2, b2, dn2)
                            if self.iflg:
                                self.matrix_S[m, n] += self.F1(a1, b1, dn2)
                                b1 -= dn
                        else :
                            # CALCULATION WITH DIELECTRICS
                            # Imn= sinm*(Ix-Ix1)-cosm*(Iy-Iy1)
                            f2= self.F2(a2, b2, dn2)
                            f3= self.F3(a2, b2, dn2)
                            Imn= ((dx   -b2*cosn)*f2-cosn*f3)*sinm -((ym-yn-b2*sinn)*f2-sinn*f3)*cosm 
                            if self.iflg:
                                f2= self.F2(a1, b1, dn2)
                                f3= self.F3(a1, b1, dn2)
                                Imn += -((dx   -b1*cosn)*f2-cosn*f3)*sinm + ((ym+yn+b1*sinn)*f2-sinn*f3)*cosm 
                                b1-= dn
                            if m==n:
                                self.diag_S11_C[n-self.n_c]=Imn + er_plus
                                self.diag_S11_L[n-self.n_c]=Imn + mu_plus
                            # HACK: This code is meaningless
                            #if Imn == 0.0:
                            #    self.matrix_S[m, n] = -1.0
                            else:
                                self.matrix_S[m, n]= Imn
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
        # FIXME: fill additional row an column for calculate L
        if not self.iflg : # fill in additional row and column
            sz = m_size-1
            n = 0
            for bound_m in self.list_bounds: # cycle through conductors
                section_m=bound_m['section']
                n_subint_m=bound_m['n_subint']
                for si in xrange(n_subint_m):
                    if bound_m['mat_type']==False:
                        si_len = section_m.len/n_subint_m
                        self.matrix_S[n, sz]=si_len/self.matrix_S[n, n]
                        self.matrix_S[sz, n]=si_len*bound_m['mat_param'].get('erp', 1.0)
                    else: # clear rest of matrix cells
                        self.matrix_S[n,sz],self.matrix_S[sz,n] = 0.0,0.0
                    n+=1

    # TODO: Refactoring
    # FIXME: diag_S11_L and diag_S11_C
    def SmnOrtho(self):
        self.matrix_S=numpy.zeros((self.m_size,self.m_size))
        bUpdate=False
        m = 0
        for bound_m in self.list_bounds: # BEGIN cycle through intervals
            section_m=bound_m['section']
            n_subint_m=bound_m['n_subint']
            # check if conductors was processed
            bDiel = bound_m['mat_type']
            # set parameters equal for all subsections
            sinm = section_m.sint
            negcosm = -section_m.cost 
            dxm,dym=section_m.dx/n_subint_m, section_m.dy/n_subint_m
            first_subsection=section_m.getSubinterval(0,n_subint_m)
            # xm, ym - centers of subinterval
            xm,ym = first_subsection.center.x, first_subsection.center.y
            er_plus=0.0
            erp = bound_m['mat_param'].get('erp', 1.0)
            erm = bound_m['mat_param'].get('erm', 1.0)
            if erp==erm:
                raise ValueError
            er_plus=(erp+erm)*pi/(erp-erm)
            # BEGIN cycle through subintervals
            for subint_m in xrange(n_subint_m):
                # DO JUST THE SAME CALCULATIONS FOR INTEGRAL INTERVALS
                n = 0
                for bound_n in self.list_bounds:
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
                                fi= -self.FI(a1, a2, c1, dn)
                                if self.iflg :
                                    if sinn==0.0 :
                                        fi+=self.FI(a1, a2, c2, dn)
                                    else :
                                        fi+=self.FI(a1+2*ym, a2+2*ym, c2, dn)
                            else :
                                # CALCULATION WITH DIELECTRICS
                                if sinm==0.0: #  m- ortho Y
                                    if sinn==0.0:
                                        # both ortho Y
                                        if c1:
                                            fi= self.sumatan(a1,-a2,c1)
                                        else:
                                            fi= 0.0
                                        if self.iflg and c2!=0.0:
                                            fi-= self.sumatan(a1,-a2,c2)
                                    
                                    else :  # m- ortho Y, n- ortho X
                                        tmp3= c1*c1; #HERE с1==с2
                                        fi= 0.5*log((a2*a2+tmp3)/(a1*a1+tmp3))
                                        if self.iflg :
                                            tmp1= a1+2*ym
                                            tmp2= tmp1- dn
                                            fi-= 0.5*log((tmp2*tmp2+tmp3)/(tmp1*tmp1+tmp3))
                                    # HACK: subintervals direction fix by T.R. Gazizov - DISABLED
                                    fi*= negcosm
                                else : #  m- ortho X
                                    if sinn==0.0 : # m- ortho X, n- ortho Y
                                        fi= 0.5*log((a2*a2+c1*c1)/(a1*a1+c1*c1))
                                        if self.iflg :
                                            fi-= 0.5*log((a2*a2+c2*c2)/(a1*a1+c2*c2))
    
                                    else :  # both ortho X
                                        if c1 :
                                            fi= self.sumatan(a1,-a2,c1)
                                        else :
                                            fi= 0.0
                                        if self.iflg and c2!=0.0:
                                            fi-= self.sumatan(a1+2*ym,-a2-2*ym,c2)
                                    # HACK: subintervals direction fix by T.R. Gazizov - DISABLED
                                    fi*= sinm
                                if(m==n):
                                    fi+= er_plus
                             # END DIELECTRIC PROCESSING
                            
                            self.matrix_S[m, n] = fi
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
        if not self.iflg : # fill in additional row and column
            sz = m_size-1
            n = 0
            for bound_m in self.list_bounds: # cycle through conductors
                section_m=bound_m['section']
                n_subint_m=bound_m['n_subint']
                for si in xrange(n_subint_m):
                    if bound_m['mat_type']==False:
                        si_len = section_m.len/n_subint_m
                        self.matrix_S[n, sz]=si_len/self.matrix_S[n, n]
                        self.matrix_S[sz, n]=si_len*bound_m['mat_param'].get('erp', 1.0)
                    else: # clear rest of matrix cells
                        self.matrix_S[n,sz],self.matrix_S[sz,n] = 0.0,0.0
                    n+=1
class LRCG(Smn):
    def calcLC(self):
        cond_sect=(filter(lambda x: x['mat_type']==False,self.list_bounds))
        n_cond=len(set(map(lambda x: x['mat_count'],cond_sect)))
        self.SmnAny2D()

        # Excitation vector filling
        exc_v = numpy.zeros((self.m_size,n_cond))
        beg,n,old_cond=0,0,cond_sect[0]['mat_count']
        for bound in cond_sect:
            if old_cond!=bound['mat_count']: n+=1
            old_cond=bound['mat_count']
            end=beg+bound['n_subint']
            exc_v[beg:end,n]=Coef_C
            beg=end

        # Matrix Q calculating
        for i in xrange(self.n_c,self.m_size):
            self.matrix_S[i,i]=self.diag_S11_C[i-self.n_c]
        self.matrix_QC=la.solve(self.matrix_S,exc_v)

        for i in xrange(self.n_c,self.m_size):
            self.matrix_S[i,i]=self.diag_S11_L[i-self.n_c]
        self.matrix_QL=la.solve(self.matrix_S,exc_v)

        # Matrix C and L calculating
        self.mC = numpy.zeros((n_cond,n_cond))
        self.mL = numpy.zeros((n_cond,n_cond))
        beg,m,old_cond=0,0,cond_sect[0]['mat_count']
        for bound in cond_sect:
            if old_cond!=bound['mat_count']: m+=1
            old_cond=bound['mat_count']
            end=beg+bound['n_subint']
            # FIXME: It works only for equable segmentation
            # For smart segmentation is necessary to take length of the every subinterval
            subint_len=bound['section'].getSubinterval(n=bound['n_subint']).len
            self.matrix_QC[beg:end,0:n_cond]*=subint_len*bound['mat_param']['erp']
            self.matrix_QL[beg:end,0:n_cond]*=subint_len/bound['mat_param']['mup']
            for n in xrange(n_cond):
                self.mC[m,n]+=self.matrix_QC[beg:end,n].sum()
                self.mL[m,n]+=self.matrix_QL[beg:end,n].sum()
            beg=end
        self.mL=la.inv(self.mL)/(V0*V0)

        for i in xrange(self.n_c,self.m_size):
            self.matrix_S[i,i]=0.0

