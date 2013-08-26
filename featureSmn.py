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
def Smn(config):
    if type(config) is Conf:
        for bound_m in config:
            erp = bound_m['mat_param'].get('erp', 1.0)
            erm = bound_m['mat_param'].get('erm', 1.0)
            tdp = bound_m['mat_param'].get('tdp', 0.0)
            tdm = bound_m['mat_param'].get('tdm', 0.0)
        #return []
    else: raise TypeError
