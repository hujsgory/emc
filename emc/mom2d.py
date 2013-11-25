#!/usr/bin/python
#coding:  utf8
from math import *
import numpy
import numpy.linalg as la
import _smn

eps0 = 8.854187817e-12 # dielectric constant
Coef_C = 4*pi*eps0
V0 = 299792458.0     # light velocity


## \class Coord
class Coord(object):

    ## \fn __init__
    # \brief Constructor, set coordinates x and y
    # \param x
    # \param y
    def __init__(self, x=0.0, y=0.0):
        self.x, self.y = float(x), float(y)

    def __eq__(self, _coord):
        return self.x == _coord.x and self.y == _coord.y

    def __ne__(self, _coord):
        return self.x != _coord.x or self.y != _coord.y

    def __str__(self):
        return str((self.x, self.y))

    def sint(self, _coord):
        return (_coord.y-self.y)/hypot(_coord.y-self.y, _coord.x-self.x)

    def cost(self, _coord):
        return (_coord.x-self.x)/hypot(_coord.y-self.y, _coord.x-self.x)

    def center(self, _coord):
        return Coord((_coord.x + self.x)/2.0, (_coord.y + self.y)/2.0)

    def len(self, _coord):
        return hypot(_coord.y-self.y, _coord.x-self.x)


## \class Section
class Section(object):

    ## \fn __init__
    # \brief Constructor, set points beg and end
    # \param beg Coord object, begin of the section
    # \param end Coord object, end of the section
    def __init__(self, beg=Coord(), end=Coord()):
        if type(beg) is Coord and type(end) is Coord:
            self.beg, self.end = beg, end
        else:
            raise TypeError

    def __eq__(self, _section):
        return self.beg == _section.beg and self.end == _section.end

    def __ne__(self, _section):
        return self.beg != _section.beg or self.end != _section.end

    def __str__(self):
        return str(str(self.beg) + '-' + str(self.end))

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

    def getSubinterval(self, i=0, n=1):
        if type(i) is int and type(n) is int:
            if i >= n or i < 0:
                raise ValueError
            x, y = self.beg.x, self.beg.y
            dx, dy = self.dx, self.dy
            return Section(Coord(x + i*dx/n, y + i*dy/n),
                           Coord(x + (i+1)*dx/n, y + (i+1)*dy/n))
        else:
            raise TypeError


## \class Structure
# \brief Object contain a list of a conductors and a dielectrics
class Structure(object):
    def __init__(self):
        self.list_cond = list()
        self.list_diel = list()
        self.iflg = True
        self.grounded = False
        self.mat_type = False     # mat_type: False - Conductor-Dielectric bound,  True - Dielectric-Dielectric bound
        self.obj_count = 0
        self.sect_count = 0
        self.mat_param = dict()
        self.iter_idx = 0

    def __iter__(self):
        self.iter_idx = 0
        return self

    def next(self):
        self.iter_idx += 1
        if self.iter_idx > len(self.list_cond) + len(self.list_diel):
            raise StopIteration
        if self.iter_idx <= len(self.list_cond):
            return self.list_cond[self.iter_idx-1]
        else:
            return self.list_diel[self.iter_idx-len(self.list_cond)-1]

    ## \fn intersection
    # \brief The function checking the two boundaries on intersection
    # \param sect1 \param sect2 two object of Section type
    def intersection(self, sect1, sect2):
        # Coeficients of Ax+By+D=0
        a1 = -sect1.dy
        b1 =  sect1.dx
        d1 = -(a1*sect1.beg.x + b1*sect1.beg.y)
        a2 = -sect2.dy
        b2 =  sect2.dx
        d2 = -(a2*sect2.beg.x + b2*sect2.beg.y)

        # Calculate halfplane of sections ends
        seg1_line2_start = a2*sect1.beg.x + b2*sect1.beg.y + d2
        seg1_line2_end   = a2*sect1.end.x + b2*sect1.end.y + d2
        seg2_line1_start = a1*sect2.beg.x + b1*sect2.beg.y + d1
        seg2_line1_end   = a1*sect2.end.x + b1*sect2.end.y + d1

        h1 = seg1_line2_start*seg1_line2_end
        h2 = seg2_line1_start*seg2_line1_end

        # Ends located in the different halfplane
        if h1 < 0. and h2 < 0.:
            return True
        # Collinear (both ends of segments lie on one line)
        if atan2(a1, b1) == atan2(a2, b2) and d1 == d2:
            #   |------| ==  ==  == |------|
            # fmin1  fmin2  fmax1  fmax2
            fmin1 = min(sect1.beg.x, sect1.end.x)
            fmin2 = min(sect2.beg.x, sect2.end.x)
            fmax1 = max(sect1.beg.x, sect1.end.x)
            fmax2 = max(sect2.beg.x, sect2.end.x)
            if fmin1 < fmax2 and fmin2 < fmax1:
                return True
            fmin1 = min(sect1.beg.y, sect1.end.y)
            fmin2 = min(sect2.beg.y, sect2.end.y)
            fmax1 = max(sect1.beg.y, sect1.end.y)
            fmax2 = max(sect2.beg.y, sect2.end.y)
            if fmin1 < fmax2 and fmin2 < fmax1:
                return True
        return False
    
    ## \fn check_intersection
    # \brief The function checking the added boundary on intersection with other boundaries 
    def check_intersection(self, section):
        return reduce(lambda r, bound:
                      r or self.intersection(bound['section'], section),
                      self.list_cond, False) or \
               reduce(lambda r, bound:
                      r or self.intersection(bound['section'], section),
                      self.list_diel, False)

    ## \fn add
    # \brief Boundary addition
    # \param section - Section() object
    # \param n_subint - number of subintervals
    def add(self, section, n_subint=1):
        if type(section) is Section and type(n_subint) is int:
            if self.check_intersection(section):
                raise ValueError
            bound = {'section': section,
                     'n_subint': n_subint,
                     'mat_param': self.mat_param,
                     'obj_count': self.obj_count,
                     'sect_count': self.sect_count}
            if self.mat_type:
                self.list_diel.append(bound)
            else:
                bound['grounded'] = self.grounded
                self.list_cond.append(bound)
            self.sect_count += 1
        else:
            raise TypeError

    ## \fn cond
    # \brief Conductor properties setting
    # \param mat_param - dictionary which can have following entries:
    #  erp - relative permittivity on right side of section,  erm - on left side
    #  tdp, tdm - dielectric loss tangent
    #  mup, mum - relative magnetic permeability
    def cond(self, **mat_param):
        self.grounded = False
        self.mat_param = mat_param
        self.mat_type = False
        self.obj_count += 1
        self.sect_count = 0

    def ground_cond(self, **mat_param):
        self.grounded = True
        self.mat_param = mat_param
        self.mat_type = False
        self.obj_count += 1
        self.sect_count = 0

    def diel(self, **mat_param):
        self.mat_param = mat_param
        self.mat_type = True
        self.obj_count += 1
        self.sect_count = 0

    ## \function set_subintervals
    #  \brief set number of a segments for the each boundary of the structure
    #  \param n_subint number of the segments
    def set_subintervals(self, n_subint):
        for bound in self:
            bound['n_subint'] = n_subint

    ## \function autosegment
    #  \brief calculate and set number of a segments to length subinterval for the each boundary of the structure
    #  \param 
    def auto_segment(self, length):
        pass

    def adaptive_segment(self, criterion):
        pass

    def smart_segment(self):
        pass

class Board(object):
    def __init__(self):
        self.layers = list()
        self.medium = {'er': 1.0, 'td': 0.0, 'mu': 1.0}

    def layer(self, height, er, td=0.0, mu=1.0):
        if height <= 0.0 or er < 1.0 or td <0.0:
            raise ValueError
        if len(self.layers) > 0:
            if len(self.layers[-1]['cond']) > 0 and \
               height <= max(map(lambda x:
                                 x['thickness']-x['depth'],
                                 self.layers[-1]['cond'])):
                raise ValueError('Thickness of conductor of previous layer is greater than height')
            if self.layers[-1]['is_cover']:
                raise ValueError('Layer can\'t be applied to cover')
        self.layers.append({'height': height,
                            'er': er,
                            'td': td,
                            'mu': mu,
                            'is_cover': False,
                            'cond': list()})

    def conductor(self, space, width, thickness, depth=0.0, grounded=False, to_center=False):
        if len(self.layers) == 0:
            raise ValueError('It is necessary to create at least one layer')
        last_layer_cond = self.layers[-1]['cond']
        if to_center:
            space -= width/2.0
            if len(last_layer_cond) > 0:
                space -= last_layer_cond[-1]['width']/2.0
        if space <= 0. or width <= 0. or thickness <= 0. or depth < 0.:
            raise ValueError('All parameters must be positive')
        if depth >= self.layers[-1]['height']:
            raise ValueError('Depth is greater than layer height')
        last_layer_cond.append({'space': space,
                                'width': width,
                                'thickness': thickness,
                                'grounded': grounded,
                                'depth': depth})

    def cover(self, height, er, td=0.0, mu=1.0):
        if len(self.layers) == 0:
            raise ValueError('It is necessary to create at least one layer')
        self.layers.append({'height': height,
                            'er': er,
                            'td': td,
                            'mu': mu,
                            'is_cover': True,
                            'cond': list()})

# FIXME: need to add approximate comparison for float numbers
    def to_structure(self):
        structure = Structure()
        # Calculation of structure's right coordinate x
# FIXME: x['cond'][0]['space'] may be raised exception
        max_x = max(map(lambda layer:
                           reduce(lambda r, cond:
                                  r + cond['space'] + cond['width'],
                                  layer['cond'],
                                  layer['cond'][0]['space']
                                  ),
                      filter(lambda layer:
                             len(layer['cond']) > 0,
                             self.layers
                             )
                      )
                  )
        # Layers drawing
        y_layer = 0.0
        for i, layer in enumerate(self.layers):
            diel_sect = list()
            y_layer += layer['height']
            er_top = self.medium['er']
            td_top = self.medium['td']
            mu_top = self.medium['mu']
            if i+1 < len(self.layers):
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
                    if not cond['grounded']:
                        structure.cond(erp=er_bottom, tdp=td_bottom, mup=mu_bottom)
                    else:
                        structure.ground_cond(erp=er_bottom, tdp=td_bottom, mup=mu_bottom)
                    y_cond_bottom = y_layer       - cond['depth']
                    y_cond_top    = y_cond_bottom + cond['thickness']
                    x_cond_left  +=                 cond['space']
                    x_cond_right  = x_cond_left   + cond['width']
                    if y_cond_bottom < y_layer:
                        beg = None
                        if y_cond_top < y_layer: beg = Coord(x_cond_left,   y_cond_top)
                        else: beg = Coord(x_cond_left, y_layer)
                        end = Coord(x_cond_left, y_cond_bottom)
                        structure.add(Section(beg, end))
                        beg = Coord(x_cond_right, y_cond_bottom)
                        if y_cond_top < y_layer:
                            end = Coord(x_cond_right, y_cond_top)
                        else:
                            end = Coord(x_cond_right, y_layer)
                        structure.add(Section(beg, end))
                    beg = Coord(x_cond_left,   y_cond_bottom)
                    end = Coord(x_cond_right,  y_cond_bottom)
                    structure.add(Section(beg, end))
                    if y_cond_top >= y_layer:
                        structure.mat_param = {'erp': er_top, 'tdp': td_top, 'mup': mu_top}
                    if y_cond_top >  y_layer:
                        beg = Coord(x_cond_right,  y_layer)
                        end = Coord(x_cond_right,  y_cond_top)
                        structure.add(Section(beg, end))
                        beg = Coord(x_cond_left,   y_cond_top)
                        end = Coord(x_cond_left,   y_layer)
                        structure.add(Section(beg, end))
                    beg = Coord(x_cond_right,  y_cond_top)
                    end = Coord(x_cond_left,   y_cond_top)
                    structure.add(Section(beg, end))
                    if cond['depth'] > cond['thickness']:
                        diel_sect.append((x_cond_left - cond['space'],  x_cond_right))
                    else:
                        diel_sect.append((x_cond_left - cond['space'],  x_cond_left))
                    x_cond_left = x_cond_right
                # Dielectric-dielectric bounds building
                structure.diel(erp=er_bottom, tdp=td_bottom, mup=mu_bottom, erm=er_top, tdm=td_top, mum=mu_top)
                for sect in diel_sect:
                    beg = Coord(sect[0],  y_layer)
                    end = Coord(sect[1],  y_layer)
                    structure.add(Section(beg, end))
                beg = Coord(x_cond_right, y_layer)
                end = Coord(max_x, y_layer)
                structure.add(Section(beg, end))
            # Cover building
            else:
                layer['cover'] = list()
                width_all = 0.0
                if self.layers[i-1]['is_cover']:
                    surface = self.layers[i-1]['cover']
                else:
                    surface = self.layers[i-1]['cond']
                for section in surface:
                    space = section.get('space', 0.0)
                    if space > 0.0:
                        layer['cover'].append({'width': space, 'thickness': 0.0})
                        width_all += section['space']
                    thickness = section['thickness']-section.get('depth', 0.0)
                    if thickness <= 0.0: thickness = 0.0
                    width_all += section['width']
                    layer['cover'].append({'width': section['width'], 'thickness': thickness})
                layer['cover'].append({'width': max_x-width_all, 'thickness': 0.0})

                cover = layer['cover']
                for j in xrange(len(cover)-1):
                    if cover[j]['thickness'] > cover[j+1]['thickness']:
                        cover[j]['width'] += layer['height']
                        cover[j+1]['width'] -= layer['height']
                    elif cover[j]['thickness'] < cover[j+1]['thickness']:
                        cover[j]['width'] -= layer['height']
                        cover[j+1]['width'] += layer['height']

                # Remove sections,  when width is negative or 0
                check_width = True
                while check_width:
                    check_width = False
                    j = 0
                    while j < len(cover):
                        if cover[j]['width'] <= 0. and len(cover) > 1:
                            if j == 0:
                                cover[j+1]['width'] += cover[j]['width']
                            elif j == len(cover)-1:
                                cover[j-1]['width'] += cover[j]['width']
                            else:
                                cover[j+1]['width'] += cover[j]['width']/2.0
                                cover[j-1]['width'] += cover[j]['width']/2.0
                            del cover[j]
                            check_width = True
                            continue
                        j += 1

                structure.diel(erp=er_bottom,
                                    tdp=td_bottom,
                                    mup=mu_bottom,
                                    erm=er_top,
                                    tdm=td_top,
                                    mum=mu_top)
                x_left,x_right = 0.0,0.0
                for section_cur, section_next in zip(cover, cover[1:]):
                    x_right = x_left + section_cur['width']
                    beg = Coord(x_left,  y_layer + section_cur['thickness'])
                    end = Coord(x_right, y_layer + section_cur['thickness'])
                    structure.add(Section(beg, end))
                    if section_cur['thickness'] != section_next['thickness']:
                        beg = Coord(x_right, y_layer + section_cur ['thickness'])
                        end = Coord(x_right, y_layer + section_next['thickness'])
                        structure.add(Section(beg, end))
                    x_left = x_right
                x_right = x_left + cover[-1]['width']
                beg = Coord(x_left,  y_layer + cover[-1]['thickness'])
                end = Coord(x_right, y_layer + cover[-1]['thickness'])
                structure.add(Section(beg, end))
        return structure


## \class Matrix
#  \brief Container for blocks of matrix S and related methods
class Matrix(object):
    def __init__(self, A00, A01, A10, A11):
        if type(A00) is numpy.ndarray:
            self.A00 = A00
            self.nc = A00.shape[0]
        else:
            raise ValueError
        if type(A01) is numpy.ndarray and type(A10) is numpy.ndarray and type(A11) is numpy.ndarray:
            if A00.shape[0] == A01.shape[0] and A10.shape[0] == A11.shape[0] and A00.shape[1] == A10.shape[1] and A01.shape[1] == A11.shape[1]:
                self.A01, self.A10, self.A11 = A01, A10, A11
                self.nd=A01.shape[1]
            else:
                self.nd = 0

    def factorize(self):
        if self.nd > 0:
            self.A10 = numpy.dot(self.A10, self.A00)
            self.A11 -= numpy.dot(self.A10, self.A01)
            self.A11 = la.inv(self.A11)

    def solve(self, b):
        nc = self.nc
        nd = self.nd
        x = b.copy()
        if self.nd > 0:
            x[nc:] -= numpy.dot(self.A10, x[:nc])
            x[nc:]  = numpy.dot(self.A11, x[nc:])
            x[:nc] -= numpy.dot(self.A01, x[nc:])
        x[:nc] = numpy.dot(self.A00, x[:nc])
        return x

    # c = A*b
    def dot(self, b, c=None):
        nc = self.nc
        nd = self.nd
        if c == None:
            c = numpy.zeros((nc + nd, b.shape[1]))
        if nd > 0:
            c[:nc] = numpy.dot(self.A00, b[:nc]) + numpy.dot(self.A01, b[nc:])
            c[nc:] = numpy.dot(self.A10, b[:nc]) + numpy.dot(self.A11, b[nc:])
        else:
            c[:] = numpy.dot(self.A00, b)
        return c

    ## \function iterative
    #  \brief Stabilized bi-conjugate gradient method with preconditioning (BiCGStab)
    #  \param M Matrix object with factorized matrixes S
    #  \param b vector of right-hand members
    def iterative(self, M, b, tol = 1e-16, max_iter = 30):
        nc = self.nc
        nd = self.nd
        n_cond = b.shape[1]
        A = self
        X = numpy.ones((nc+nd, n_cond))
        V = numpy.zeros((nc+nd, n_cond))
        P = numpy.zeros((nc+nd, n_cond))
        R = b - A.dot(X)
        Rt = R.copy()
        S = numpy.zeros((nc+nd, n_cond))
        T = numpy.zeros((nc+nd, n_cond))
        normR0 = la.norm(R)
        alpha = numpy.ones(n_cond)
        beta = numpy.zeros(n_cond)
        rho = numpy.zeros(n_cond)
        rho_old = numpy.ones(n_cond)
        omega = numpy.ones(n_cond)
        for iter in xrange(max_iter):
            for i in xrange(n_cond):
                rho[i] = numpy.dot(Rt[:,i], R[:,i])
                if rho[i] == 0.0:
                    raise ValueError
                beta[i] = (rho[i]/rho_old[i])*(alpha[i]/omega[i])
            P = R + beta*(P - omega*V)
            Pt = M.solve(P)
            A.dot(Pt, V)
            for i in xrange(n_cond):
                alpha[i] = rho[i]/numpy.dot(Rt[:,i], V[:,i])
                S[:,i] = R[:,i] - alpha[i]*V[:,i]
                X[:,i] += alpha[i]*Pt[:,i]
            if la.norm(S)/normR0 <= tol:
                break
            St = M.solve(S)
            A.dot(St, T)
            for i in xrange(n_cond):
                omega[i] = numpy.dot(T[:,i], S[:,i])/numpy.dot(T[:,i], T[:,i])
                X[:,i] += omega[i]*St[:,i]
                R[:,i] = S[:,i] - omega[i]*T[:,i]
            if la.norm(R)/normR0 <= tol:
                 break
            rho_old = rho
        return X


## \class Smn
# \brief Matrix which binds a vector of charges and a vector of potential
class Smn(object):
    def __init__(self, conf):
        if type(conf) is not Structure:
            raise TypeError
        self.list_cond = conf.list_cond
        self.list_diel = conf.list_diel
        self.not_grounded_cond = filter(lambda x: not x['grounded'], self.list_cond)
        if len(self.not_grounded_cond) <= 0:
            raise ValueError('Not grounded conductors is not exist')
        self.iflg = conf.iflg
        self.nc = reduce(lambda r, x: r + x['n_subint'], self.list_cond, 0)
        self.n_cond = len(set(map(lambda x: x['obj_count'], self.not_grounded_cond)))
        self.exc_v0 = numpy.zeros((self.nc, self.n_cond))
        beg, n, old_cond = 0, 0, self.not_grounded_cond[0]['obj_count']
        for bound in self.list_cond:
            end = beg + bound['n_subint']
            if not bound['grounded']:
                if old_cond != bound['obj_count']:
                    n += 1
                old_cond = bound['obj_count']
                self.exc_v0[beg: end, n] = Coef_C
            beg = end
        self.S00 = numpy.zeros((self.nc, self.nc), dtype=numpy.float64)
        self._calcSmn_(self.S00, self.list_cond, self.list_cond, False)
        self.diag_S00 = self.S00.diagonal()
        self.is_inv_S00 = False

    def _calcSmn_(self, block_S, list1, list2, bDiel):
        for bounds in self.list_cond:
            section = bounds['section']
            bounds['_section_'] = [section.beg.x, section.beg.y, section.end.x, section.end.y]
        for bounds in self.list_diel:
            section = bounds['section']
            bounds['_section_'] = [section.beg.x, section.beg.y, section.end.x, section.end.y]
        if type(block_S) is not numpy.ndarray or \
           not numpy.issubdtype(block_S.dtype, numpy.float64) or \
           type(list1) is not list or \
           type(list2) is not list:
            raise TypeError
        if len(block_S.shape) == 2 and \
           block_S.shape[0] == reduce(lambda r, x: r + x['n_subint'], list1, 0) and \
           block_S.shape[1] == reduce(lambda r, x: r + x['n_subint'], list2, 0):
            _smn._calcSmn_(block_S, list1, list2, bDiel, self.iflg)
        else:
            raise ValueError

    @property
    def list_diel_C(self):
        return filter(lambda x:
                      x['mat_param'].get('erp', 1.0) != x['mat_param'].get('erm', 1.0),
                      self.list_diel)

    @property
    def list_diel_L(self):
        return filter(lambda x:
                      x['mat_param'].get('mup', 1.0) != x['mat_param'].get('mum', 1.0),
                      self.list_diel)

    @property
    def nd_C(self):
        nd = 0
        if not self.iflg:
            nd = 1
        return reduce(lambda r, x: r + x['n_subint'], self.list_diel_C, nd)

    @property
    def nd_L(self):
        nd = 0
        if not self.iflg:
            nd = 1
        return reduce(lambda r, x: r + x['n_subint'], self.list_diel_L, nd)

    def fill_SC(self):
        nc=self.nc
        nd=self.nd_C
        S01_C = numpy.zeros((nc, nd), dtype=numpy.float64)
        S10_C = numpy.zeros((nd, nc), dtype=numpy.float64)
        S11_C = numpy.zeros((nd, nd), dtype=numpy.float64)
        list_cond = self.list_cond
        list_diel = self.list_diel_C
        self._calcSmn_(S01_C, list_cond, list_diel, False)
        self._calcSmn_(S10_C, list_diel, list_cond, True)
        self._calcSmn_(S11_C, list_diel, list_diel, True)
        m = 0
        for bound in list_diel:
            erp = bound['mat_param'].get('erp', 1.0)
            erm = bound['mat_param'].get('erm', 1.0)
            er = (erp + erm) * pi / (erp - erm)
            for i in xrange(bound['n_subint']):
                S11_C[m, m] += er
                m += 1
        self.matrix_SC = Matrix(self.S00, S01_C, S10_C, S11_C)
        if not self.iflg:
            lastC = self.nd_C-1
            n = 0
            for bound in self.list_cond:
                section_m = bound['section']
                n_subint_m = bound['n_subint']
                for si in xrange(n_subint_m):
                    subint_len = bound['section'].getSubinterval(i, n_subint_m).len
                    self.matrix_SC.A01[n, lastC] = subint_len/self.diag_S00[n]
                    self.matrix_SC.A10[lastC, n] = subint_len*bound_m['mat_param'].get('erp',  1.0)
                    n += 1
    
    def fill_SL(self):
        nc=self.nc
        nd=self.nd_L
        S01_L = numpy.zeros((nc, nd), dtype=numpy.float64)
        S10_L = numpy.zeros((nd, nc), dtype=numpy.float64)
        S11_L = numpy.zeros((nd, nd), dtype=numpy.float64)
        list_cond = self.list_cond
        list_diel = self.list_diel_L
        self._calcSmn_(S01_L, list_cond, list_diel, False)
        self._calcSmn_(S10_L, list_diel, list_cond, True)
        self._calcSmn_(S11_L, list_diel, list_diel, True)
        m = 0
        for bound in list_diel:
            mup = bound['mat_param'].get('mup', 1.0)
            mum = bound['mat_param'].get('mum', 1.0)
            mu = (mup + mum) * pi / (mum - mup)
            for i in xrange(bound['n_subint']):
                S11_L[m, m] += mu
                m += 1
        self.matrix_SL = Matrix(self.S00, S01_L, S10_L,  S11_L)
        if not self.iflg:
            lastL = self.nd_L-1
            n = 0
            for bound in self.list_cond:
                section_m = bound['section']
                n_subint_m = bound['n_subint']
                for si in xrange(n_subint_m):
                    subint_len = bound['section'].getSubinterval(i, n_subint_m).len
                    self.matrix_SL.A01[n, lastL] = subint_len/self.diag_S00[n]
                    self.matrix_SL.A10[lastL, n] = subint_len/bound_m['mat_param'].get('mup',  1.0)
                    n += 1

    def factorize_SC(self):
        if not self.is_inv_S00:
            self.S00[:] = la.inv(self.S00)
            self.is_inv_S00 = True
        self.matrix_SC.factorize()
    
    def factorize_SL(self):
        if not self.is_inv_S00:
            self.S00[:] = la.inv(self.S00)
            self.is_inv_S00 = True
        self.matrix_SL.factorize()

    def solve_SC(self):
        matrix_Q = numpy.zeros((self.nc + self.nd_C, self.n_cond))
        matrix_Q[:self.nc] = self.exc_v0
        return self.matrix_SC.solve(matrix_Q)

    def solve_SL(self):
        matrix_Q = numpy.zeros((self.nc + self.nd_L, self.n_cond))
        matrix_Q[:self.nc] = self.exc_v0
        return self.matrix_SL.solve(matrix_Q)

    def iterative_C(self, M):
        matrix_Q = numpy.zeros((self.nc + self.nd_C, self.n_cond))
        matrix_Q[:self.nc] = self.exc_v0
        return self.matrix_SC.iterative(M.matrix_SC, matrix_Q)

    def iterative_L(self, M):
        matrix_Q = numpy.zeros((self.nc + self.nd_L, self.n_cond))
        matrix_Q[:self.nc] = self.exc_v0
        return self.matrix_SL.iterative(M.matrix_SL, matrix_Q)


class RLGC(object):
    def __init__(self, conf):
        self.is_calc_C, self.is_calc_L = False, False
        self.smn = Smn(conf)

    def update(self, conf):
        self.precondition = self.smn
        self.smn = Smn(conf)
        if self.is_calc_C:
            self.smn.fill_SC()
            self.matrix_QC = self.smn.iterative_C(self.precondition)
            self.mC[:,:] = 0.0
            self.calc_C()
        if self.is_calc_L:
            self.smn.fill_SL()
            self.matrix_QL = self.smn.iterative_L(self.precondition)
            self.mL[:,:] = 0.0
            self._calc_L()
    
    def calc_C(self):
        self.is_calc_C = True
        self.smn.fill_SC()
        self.smn.factorize_SC()
        self.matrix_QC = self.smn.solve_SC()
        self.mC = numpy.zeros((self.smn.n_cond, self.smn.n_cond))
        self._calc_C()
    
    def _calc_C(self):
        beg, m, old_cond = 0, 0, self.smn.not_grounded_cond[0]['obj_count']
        for bound in self.smn.list_cond:
            end = beg + bound['n_subint']
            if not bound['grounded']:
                if old_cond != bound['obj_count']:
                    m += 1
                old_cond = bound['obj_count']
            erp = bound['mat_param'].get('erp', 1.0)
            if not bound['grounded']:
                for j, i in enumerate(xrange(beg, end)):
                    subint_len = bound['section'].getSubinterval(j, bound['n_subint']).len
                    self.matrix_QC[i, 0:self.smn.n_cond] *= subint_len*erp
                for n in xrange(self.smn.n_cond):
                    self.mC[m, n] += self.matrix_QC[beg: end, n].sum()
            beg = end
    
    def calc_L(self):
        self.is_calc_L = True
        self.smn.fill_SL()
        self.smn.factorize_SL()
        self.matrix_QL = self.smn.solve_SL()
        self.mL = numpy.zeros((self.smn.n_cond, self.smn.n_cond))
        self._calc_L()
    
    def _calc_L(self):
        beg, m, old_cond = 0, 0, self.smn.not_grounded_cond[0]['obj_count']
        for bound in self.smn.list_cond:
            end = beg + bound['n_subint']
            if not bound['grounded']:
                if old_cond != bound['obj_count']:
                    m += 1
                old_cond = bound['obj_count']
            mup = bound['mat_param'].get('mup', 1.0)
            if not bound['grounded']:
                for j, i in enumerate(xrange(beg, end)):
                    subint_len = bound['section'].getSubinterval(j, bound['n_subint']).len
                    self.matrix_QL[i, 0:self.smn.n_cond] *= subint_len/mup
                for n in xrange(self.smn.n_cond):
                    self.mL[m, n] += self.matrix_QL[beg: end, n].sum()
            beg = end
        self.mL = la.inv(self.mL)/(V0*V0)

    def calc_LC(self):
        self.calc_L()
        self.calc_C()
        
