#!/usr/bin/python
#coding: utf8
from math import *
import numpy
import numpy.linalg as la

eps0=8.854187817e-12 #dielectric constant
Coef_C = 4*pi*eps0   
V0 = 299792458.0     #light velocity

class Coord(object):
    def __init__(self,x=0.0,y=0.0):
        self.x,self.y=float(x),float(y)
    def __eq__(self,_coord):
        return self.x==_coord.x and self.y==_coord.y
    def __ne__(self,_coord):
        return self.x!=_coord.x or self.y!=_coord.y
    def __str__(self):
        return str((self.x,self.y))
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
    def __str__(self):
        return str(str(self.beg)+'-'+str(self.end))
    def __repr__(self):
        return str(self)
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
        self.list_cond=list()
        self.list_diel=list()
        self.iflg=True
        self.grounded=False
        self.mat_type=False     #mat_type: False - Conductor-Dielectric bound, True - Dielectric-Dielectric bound
        self.obj_count=0
        self.sect_count=0
        self.mat_param=dict()
        self.iter_idx=0
    def __iter__(self):
        self.iter_idx=0
        return self
    def next(self):
        self.iter_idx+=1
        if self.iter_idx>len(self.list_cond)+len(self.list_diel):
            raise StopIteration
        if self.iter_idx<=len(self.list_cond):
            return self.list_cond[self.iter_idx-1]
        else:
            return self.list_diel[self.iter_idx-len(self.list_cond)-1]
    def intersection(self,sect1,sect2):
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
    def check_intersection(self,section):
        return reduce(lambda r,bound: r or self.intersection(bound['section'],section),self.list_cond,False) or \
               reduce(lambda r,bound: r or self.intersection(bound['section'],section),self.list_diel,False)
            
    # n_subint: number of bound_m's subintervals
    def add(self,section,n_subint=1):
        if type(section) is Section and type(n_subint) is int:
            if self.check_intersection(section): raise ValueError
            bound={'section':section,'n_subint':n_subint,'mat_param':self.mat_param, 'obj_count': self.obj_count,'sect_count': self.sect_count}
            if self.mat_type:
                self.list_diel.append(bound)
            else:
                bound['grounded']=self.grounded
                self.list_cond.append(bound)
            self.sect_count+=1
        else: raise TypeError
    # mat_param: erp - relative permittivity on right side of section, erm - on left side
    # tdp,tdm - dielectric loss tangent 
    # mup,mum - relative magnetic permeability
    def cond(self,**mat_param):
        self.grounded=False
        self.mat_param=mat_param
        self.mat_type=False
        self.obj_count+=1
        self.sect_count=0
    def ground_cond(self,**mat_param):
        self.grounded=True
        self.mat_param=mat_param
        self.mat_type=False
        self.obj_count+=1
        self.sect_count=0
    def diel(self,**mat_param):
        self.mat_param=mat_param
        self.mat_type=True
        self.obj_count+=1
        self.sect_count=0


class Board():
    def __init__(self):
        self.layers=list()
        self.medium={'er':1.0,'td':0.0,'mu':1.00000037}
    def layer(self,height,er,td=0.0,mu=1.0):
        if height<=0.0 or er<1.0 or td <0.0:
            raise ValueError
        if len(self.layers)>0:
            if len(self.layers[-1]['cond'])>0 and height<=max(map(lambda x: x['thickness']-x['depth'], self.layers[-1]['cond'])):
                raise ValueError('Thickness of conductor of previous layer is greater than height')
            if self.layers[-1]['is_cover']:
                raise ValueError('Layer can\'t be applied to cover')
        self.layers.append({'height':height,'er':er,'td':td,'mu':mu,'is_cover':False,'cond':list()})
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
    def cover(self,height,er,td=0.0,mu=1.0):
        if len(self.layers)==0:
            raise ValueError('It is necessary to create at least one layer')
        self.layers.append({'height':height,'er':er,'td':td,'mu':mu,'is_cover':True,'cond':list()})
# FIXME: need to add approximate comparison for float numbers
    def board2conf(self):
        self.conf=Conf()
        # Calculation of structure's right coordinate x
# FIXME: x['cond'][0]['space'] may be raised exception
        max_x=max(map(lambda layer:
                           reduce(lambda r,cond:                     \
                                      r+cond['space']+cond['width'], \
                                  layer['cond'],                     \
                                  layer['cond'][0]['space']          \
                                  ),                                 \
                      filter(lambda layer:                           \
                                 len(layer['cond'])>0,               \
                             self.layers                             \
                             )                                       \
                      )                                              \
                  )
        # Layers drawing
        y_layer=0.0
        for i,layer in enumerate(self.layers):
            diel_sect=list()
            y_layer+=layer['height']
            er_top=self.medium['er']
            td_top=self.medium['td']
            mu_top=self.medium['mu']
            if i+1<len(self.layers):
                er_top = self.layers[i+1]['er']
                td_top = self.layers[i+1]['td']
                mu_top = self.layers[i+1]['mu']
            er_bottom = layer['er']
            td_bottom = layer['td']
            mu_bottom = layer['mu']
            if not layer['is_cover']:
                x_cond_left  = 0.0
                x_cond_right = 0.0
                # Conductor-dielectric bounds building
                for cond in layer['cond']:
                    self.conf.cond(erp=er_bottom,tdp=td_bottom,mup=mu_bottom)
                    y_cond_bottom = y_layer       - cond['depth']
                    y_cond_top    = y_cond_bottom + cond['thickness']
                    x_cond_left  +=                 cond['space']
                    x_cond_right  = x_cond_left   + cond['width']
                    if y_cond_bottom < y_layer :
                        beg=None
                        if y_cond_top < y_layer : beg=Coord(x_cond_left,  y_cond_top)
                        else:                     beg=Coord(x_cond_left,  y_layer)
                        end=Coord(x_cond_left, y_cond_bottom)
                        self.conf.add(Section(beg,end))
                        beg = Coord(x_cond_right, y_cond_bottom);
                        if y_cond_top < y_layer : end=Coord(x_cond_right, y_cond_top)
                        else:                     end=Coord(x_cond_right, y_layer)
                        self.conf.add(Section(beg,end))
                    beg=Coord(x_cond_left,  y_cond_bottom)
                    end=Coord(x_cond_right, y_cond_bottom)
                    self.conf.add(Section(beg,end))
                    if y_cond_top >= y_layer: self.conf.mat_param={'erp':er_top,'tdp':td_top,'mup':mu_top}
                    if y_cond_top >  y_layer:
                        beg = Coord(x_cond_right, y_layer)
                        end = Coord(x_cond_right, y_cond_top)
                        self.conf.add(Section(beg,end))
                        beg = Coord(x_cond_left,  y_cond_top)
                        end = Coord(x_cond_left,  y_layer)
                        self.conf.add(Section(beg,end))
                
                    beg = Coord(x_cond_right, y_cond_top)
                    end = Coord(x_cond_left,  y_cond_top);
                    self.conf.add(Section(beg,end))
                    if cond['depth'] > cond['thickness']:
                        diel_sect.append((x_cond_left - cond['space'], x_cond_right))
                    else:
                        diel_sect.append((x_cond_left - cond['space'], x_cond_left))
                    x_cond_left = x_cond_right
                # Dielectric-dielectric bounds building
                self.conf.diel(erp=er_bottom,tdp=td_bottom,mup=mu_bottom,erm=er_top,tdm=td_top,mum=mu_top)
                for sect in diel_sect:
                    beg = Coord(sect[0], y_layer)
                    end = Coord(sect[1], y_layer)
                    self.conf.add(Section(beg,end))
                beg = Coord(x_cond_right, y_layer)
                end = Coord(max_x,   y_layer)
                self.conf.add(Section(beg,end))
            # Cover building
            else:
                layer['cover']=list()
                width_all=0.0
                if self.layers[i-1]['is_cover']:
                    surface=self.layers[i-1]['cover']
                else:
                    surface=self.layers[i-1]['cond']
                for section in surface:
                    space=section.get('space',0.0)
                    if space>0.0:
                        layer['cover'].append({'width':space,'thickness':0.0})
                        width_all+=section['space']
                    thickness=section['thickness']-section.get('depth',0.0)
                    if thickness<=0.0: thickness=0.0
                    width_all+=section['width']
                    layer['cover'].append({'width':section['width'],'thickness':thickness})
                layer['cover'].append({'width':max_x-width_all,'thickness':0.0})

                cover=layer['cover']
                for j in xrange(len(cover)-1):
                    if cover[j]['thickness']>cover[j+1]['thickness']:
                        cover[j  ]['width']+=layer['height']
                        cover[j+1]['width']-=layer['height']
                    elif cover[j]['thickness']<cover[j+1]['thickness']:
                        cover[j  ]['width']-=layer['height']
                        cover[j+1]['width']+=layer['height']
                
                # Remove sections, when width is negative or 0
                check_width=True
                while check_width :
                    check_width=False
                    j=0
                    while j<len(cover) :
                        if cover[j]['width']<=0. and len(cover)>1 :
                            if j==0 :
                                cover[j+1]['width']+=cover[j]['width']
                            elif j==len(cover)-1 :
                                cover[j-1]['width']+=cover[j]['width']
                            else :
                                cover[j+1]['width']+=cover[j]['width']/2.0
                                cover[j-1]['width']+=cover[j]['width']/2.0
                            del cover[j]
                            check_width=True
                            continue
                        j+=1
                # 
                self.conf.diel(erp=er_bottom,tdp=td_bottom,mup=mu_bottom,erm=er_top,tdm=td_top,mum=mu_top)
                x_left,x_right=0.0,0.0
                for section_cur,section_next in zip(cover,cover[1:]) :
                    x_right = x_left + section_cur['width']
                    beg=Coord(x_left, y_layer+section_cur['thickness'])
                    end=Coord(x_right,y_layer+section_cur['thickness'])
                    self.conf.add(Section(beg,end))
                    if section_cur['thickness']!=section_next['thickness']:
                        beg=Coord(x_right,y_layer+section_cur ['thickness'])
                        end=Coord(x_right,y_layer+section_next['thickness'])
                        self.conf.add(Section(beg,end))
                    x_left=x_right
                x_right = x_left + cover[-1]['width']
                beg=Coord(x_left, y_layer+cover[-1]['thickness'])
                end=Coord(x_right,y_layer+cover[-1]['thickness'])
                self.conf.add(Section(beg,end))

'''
Port from smn.cpp
'''

class Smn(object):
    # list_bounds,iflg,matrix_S,m_size
    def __init__(self,conf):
        if type(conf) is not Conf:
            raise TypeError
        self.list_cond=conf.list_cond
        self.list_diel=conf.list_diel
        self.iflg=conf.iflg
        self.nc=reduce(lambda r,x: r+x['n_subint'],filter(lambda x: not x['grounded'],self.list_cond),0)
        self.isCalcC,self.isCalcL=False,False

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
    
    def _calcSmn_(self,block_S,list1,list2,bDiel):
        m=0
        for bound_m in list1:
            section_m = bound_m['section']
            sinm,cosm=section_m.sint,section_m.cost
            len_bound_m=bound_m['n_subint']
            for i in xrange(len_bound_m):
                subsection_i=section_m.getSubinterval(i,len_bound_m)
                xm,ym = subsection_i.center.x,subsection_i.center.y
                n=0
                for bound_n in list2:
                    section_n = bound_n['section']
                    sinn,cosn=section_n.sint,section_n.cost
                    len_bound_n=bound_n['n_subint']
                    for j in xrange(len_bound_n):
                        subsection_j=section_n.getSubinterval(j,len_bound_n)
                        dn2 = subsection_j.len/2
                        xn,yn=subsection_j.center.x,subsection_j.center.y
                        dx = xm - xn
                        a2=dx*sinn-(ym-yn)*cosn
                        b2=dx*cosn+(ym-yn)*sinn
                        a1= dx*sinn+(ym+yn)*cosn
                        b1= dx*cosn-(ym+yn)*sinn
                        if not bDiel:
                            block_S[m, n] = -self.F1(a2, b2, dn2)
                            if self.iflg:
                                block_S[m, n] += self.F1(a1, b1, dn2)
                        else:
                            f2= self.F2(a2, b2, dn2)
                            f3= self.F3(a2, b2, dn2)
                            Imn=        ((dx   -b2*cosn)*f2-cosn*f3)*sinm - ((ym-yn-b2*sinn)*f2-sinn*f3)*cosm 
                            if self.iflg:
                                f2= self.F2(a1, b1, dn2)
                                f3= self.F3(a1, b1, dn2)
                                Imn += -((dx   -b1*cosn)*f2-cosn*f3)*sinm + ((ym+yn+b1*sinn)*f2-sinn*f3)*cosm 
                            block_S[m, n]= Imn
                        n+=1
                m+=1

    # fill in additional row and column if infinite ground is not exist
    def calcLast(self):
        if not self.iflg: 
            lastC=self.nd_C-1
            lastL=self.nd_L-1
            n = 0
            for bound in self.list_cond:
                section_m=bound['section']
                n_subint_m=bound['n_subint']
                for si in xrange(n_subint_m):
                    subint_len=bound['section'].getSubinterval(i,n_subint_m).len
                    if isCalcC:
                        self.matrix_S01_C[n, lastC]=subint_len/self.matrix_S00[n, n]
                        self.matrix_S10_C[lastC, n]=si_len*bound_m['mat_param'].get('erp', 1.0)
                    if isCalcL:
                        self.matrix_S01_L[n, lastL]=subint_len/self.matrix_S00[n, n]
                        self.matrix_S10_L[lastL, n]=si_len/bound_m['mat_param'].get('mup', 1.0)
                    n+=1
    @property
    def list_diel_C(self):
        return filter(lambda x: x['mat_param'].get('erp',1.0)!=x['mat_param'].get('erm',1.0),self.list_diel)
    @property
    def list_diel_L(self):
        return filter(lambda x: x['mat_param'].get('mup',1.0)!=x['mat_param'].get('mum',1.0),self.list_diel)
    @property
    def nd_C(self):
        nd=0
        if not self.iflg: nd=1
        return reduce(lambda r,x: r+x['n_subint'],self.list_diel_C,nd)
    @property
    def nd_L(self):
        nd=0
        if not self.iflg: nd=1
        return reduce(lambda r,x: r+x['n_subint'],self.list_diel_L,nd)
    # matrix S00 filling
    def calcS00(self):
        self.matrix_S00=numpy.zeros((self.nc,self.nc))
        self._calcSmn_(self.matrix_S00,self.list_cond,self.list_cond,False)
    # matrix S01 filling
    def calcS01(self):
        if self.isCalcC:
            self.matrix_S01_C=numpy.zeros((self.nc,self.nd_C))
            self._calcSmn_(self.matrix_S01_C,self.list_cond,self.list_diel_C,False)
        if self.isCalcL:
            self.matrix_S01_L=numpy.zeros((self.nc,self.nd_L))
            self._calcSmn_(self.matrix_S01_L,self.list_cond,self.list_diel_L,False)
# FIXME: filling last column if infinite ground is not exist
        if not self.iflg:
            pass
    # matrix S10 filling
    def calcS10(self):
        if self.isCalcC:
            self.matrix_S10_C=numpy.zeros((self.nd_C,self.nc))
            self._calcSmn_(self.matrix_S10_C,self.list_diel_C,self.list_cond,True)
        if self.isCalcL:
            self.matrix_S10_L=numpy.zeros((self.nd_L,self.nc))
            self._calcSmn_(self.matrix_S10_L,self.list_diel_L,self.list_cond,True)
# FIXME: filling last row if infinite ground is not exist
        if not self.iflg:
            pass
    # matrix S11 filling
#FIXME: Need to adding er_plus and/or mu_plus to diagonal of S11 in moment, when LC-parameters is calculated
    def calcS11(self):
        if self.isCalcC:
            nd=self.nd_C
            list_diel=self.list_diel_C
            self.matrix_S11_C=numpy.zeros((nd,nd))
            self._calcSmn_(self.matrix_S11_C,list_diel,list_diel,True)
            m=0
            for bound in list_diel:
                erp=bound['mat_param'].get('erp',1.0)
                erm=bound['mat_param'].get('erm',1.0)
                er_plus=(erp+erm)*pi/(erp-erm)
                for i in xrange(bound['n_subint']):
                    self.matrix_S11_C[m,m]+=er_plus
                    m+=1
        if self.isCalcL:
            nd=self.nd_L
            list_diel=self.list_diel_L
            self.matrix_S11_L=numpy.zeros((nd,nd))
            self._calcSmn_(self.matrix_S11_L,list_diel,list_diel,True)
            m=0
            for bound in list_diel:
                mup=bound['mat_param'].get('mup',1.0)
                mum=bound['mat_param'].get('mum',1.0)
                mu_plus=(mup+mum)*pi/(mum-mup)
                for i in xrange(bound['n_subint']):
                    self.matrix_S11_L[m,m]+=mu_plus
                    m+=1
                    
    def calcS(self):
        # A
        self.calcS00()
        self.matrix_S00=la.inv(self.matrix_S00)
        # B
        self.calcS01()
        # C
        self.calcS10()
        self.calcLast()
        if self.isCalcC:
            self.matrix_S10_C=numpy.dot(self.matrix_S10_C,self.matrix_S00)
        if self.isCalcL:
            self.matrix_S10_L=numpy.dot(self.matrix_S10_L,self.matrix_S00)
        # D
        self.calcS11()
        if self.isCalcC:
            self.matrix_S11_C-=numpy.dot(self.matrix_S10_C,self.matrix_S01_C)
            self.matrix_S11_C=la.inv(self.matrix_S11_C)
        if self.isCalcL:
            self.matrix_S11_L-=numpy.dot(self.matrix_S10_L,self.matrix_S01_L)
            self.matrix_S11_L=la.inv(self.matrix_S11_L)



class RLGC():
    def __init__(self,conf):
        self.smn=Smn(conf)
        self.n_cond=len(set(map(lambda x: x['obj_count'],self.smn.list_cond)))
    def calcC(self):
        self.smn.isCalcC=True
        self._calcLC_()
    def calcL(self):
        self.smn.isCalcL=True
        self._calcLC_()
    def calcLC(self):
        self.smn.isCalcC,self.smn.isCalcL=True,True
        self._calcLC_()

    def _calcLC_(self):
        self.smn.calcS()
        
        # Excitation vector filling
        exc_v0 = numpy.zeros((self.smn.nc,self.n_cond))
        beg,n,old_cond=0,0,self.smn.list_cond[0]['obj_count']
        for bound in self.smn.list_cond:
            if old_cond!=bound['obj_count'] and not bound['grounded']: n+=1
            old_cond=bound['obj_count']
            end=beg+bound['n_subint']
            exc_v0[beg:end,n]=Coef_C
            beg=end
        
        # Matrix Q calculating
        if self.smn.isCalcC:
            self.matrix_QC = numpy.zeros((self.smn.nc+self.smn.nd_C,self.n_cond))
            self.matrix_QC[            :self.smn.nc] = exc_v0
            self.matrix_QC[ self.smn.nc:           ]-= numpy.dot(self.smn.matrix_S10_C, self.matrix_QC[            :self.smn.nc ])
            self.matrix_QC[ self.smn.nc:           ] = numpy.dot(self.smn.matrix_S11_C, self.matrix_QC[ self.smn.nc:            ])
            self.matrix_QC[            :self.smn.nc]-= numpy.dot(self.smn.matrix_S01_C, self.matrix_QC[ self.smn.nc:            ])
            self.matrix_QC[            :self.smn.nc] = numpy.dot(self.smn.matrix_S00  , self.matrix_QC[            :self.smn.nc ])
            self.mC = numpy.zeros((self.n_cond,self.n_cond))

        if self.smn.isCalcL:
            self.matrix_QL = numpy.zeros((self.smn.nc+self.smn.nd_C,self.n_cond))
            self.matrix_QL[            :self.smn.nc] = exc_v0
            self.matrix_QL[ self.smn.nc:           ]-= numpy.dot(self.smn.matrix_S10_L, self.matrix_QL[            :self.smn.nc ])
            self.matrix_QL[ self.smn.nc:           ] = numpy.dot(self.smn.matrix_S11_L, self.matrix_QL[ self.smn.nc:            ])
            self.matrix_QL[            :self.smn.nc]-= numpy.dot(self.smn.matrix_S01_L, self.matrix_QL[ self.smn.nc:            ])
            self.matrix_QL[            :self.smn.nc] = numpy.dot(self.smn.matrix_S00  , self.matrix_QL[            :self.smn.nc ])
            self.mL = numpy.zeros((self.n_cond,self.n_cond))

        # Matrix C and L calculating
        beg,m,old_cond=0,0,self.smn.list_cond[0]['obj_count']
        for bound in self.smn.list_cond:
            if old_cond!=bound['obj_count'] and not bound['grounded']: m+=1
            old_cond=bound['obj_count']
            end=beg+bound['n_subint']
            subint_len=bound['section'].getSubinterval(n=bound['n_subint']).len
            if not bound['grounded']:
                if self.smn.isCalcC:
                    self.matrix_QC[beg:end,0:self.n_cond]*=subint_len*bound['mat_param'].get('erp',1.0)
                if self.smn.isCalcL:
                    self.matrix_QL[beg:end,0:self.n_cond]*=subint_len/bound['mat_param'].get('mup',1.0)
                for n in xrange(self.n_cond):
                    if self.smn.isCalcC: 
                        self.mC[m,n]+=self.matrix_QC[beg:end,n].sum()
                    if self.smn.isCalcL:
                        self.mL[m,n]+=self.matrix_QL[beg:end,n].sum()
            beg=end
        if self.smn.isCalcL: self.mL=la.inv(self.mL)/(V0*V0)
