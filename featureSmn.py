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
