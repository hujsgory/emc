    '''
    def test_a1(self):
        self.assertAlmostEqual(a1(self.s1,self.s2)                                    , sqrt(2.0)         , 15 , (a1(self.s1,self.s2)  ,  sqrt(2.0)))
    def test_b1(self):
        self.assertAlmostEqual(b1(self.s1,self.s2)                                    , 0.0               , 15 , (b1(self.s1,self.s2)  ,  0.0))
    def test_a2(self):
        self.assertAlmostEqual(a2(self.s1,self.s2)                                    , sqrt(0.5)         , 14 , (a2(self.s1,self.s2)  ,  sqrt(0.5)))
    def test_b2(self):
        self.assertAlmostEqual(b2(self.s1,self.s2)                                    ,-sqrt(0.5)         , 14 , (b2(self.s1,self.s2)  , -sqrt(0.5)))
    def test_F1(self):
        self.assertAlmostEqual(F1(a1(self.s1,self.s2),b1(self.s1,self.s2),self.s2.len), 1.0901906025902008, 15 , F1(a1(self.s1,self.s2), b1(self.s1,self.s2),self.s2.len))
        self.assertAlmostEqual(F1(a2(self.s1,self.s2),b2(self.s1,self.s2),self.s2.len), 0.033148387615394 , 15 , F1(a2(self.s1,self.s2), b2(self.s1,self.s2),self.s2.len))
    def test_F2(self):
        self.assertAlmostEqual(F2(a1(self.s1,self.s2),b1(self.s1,self.s2),self.s2.len), 0.655696736810798 , 15 , F2(a1(self.s1,self.s2), b1(self.s1,self.s2),self.s2.len))
        self.assertAlmostEqual(F2(a2(self.s1,self.s2),b2(self.s1,self.s2),self.s2.len), 1.565744732268386 , 15 , F2(a2(self.s1,self.s2), b2(self.s1,self.s2),self.s2.len))
        self.assertAlmostEqual(F2(                0.0,b1(self.s1,self.s2),self.s2.len), 2.0*sqrt(2)       , 15 , F2(0.0                , b1(self.s1,self.s2),self.s2.len))
        self.assertAlmostEqual(F2(                0.0,b2(self.s1,self.s2),self.s2.len), sqrt(2)           , 15 , F2(0.0                , b2(self.s1,self.s2),self.s2.len))
    def test_F3(self):
        self.assertAlmostEqual(F3(a1(self.s1,self.s2),b1(self.s1,self.s2),self.s2.len), 0.0               , 15 , F3(a1(self.s1,self.s2), b1(self.s1,self.s2),self.s2.len))
        self.assertAlmostEqual(F3(a2(self.s1,self.s2),b2(self.s1,self.s2),self.s2.len), 0.5*log(5)        , 15 , F3(a2(self.s1,self.s2), b2(self.s1,self.s2),self.s2.len))
    def test_Imn(self):
        self.assertAlmostEqual(Imn (self.s1,self.s2)                                  ,-0.213850135243756 , 15 , Imn(self.s1,self.s2))
        self.assertAlmostEqual(I_mn(self.s1,self.s2)                                  , 0.655696736810798 , 15 , I_mn(self.s1,self.s2))
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
