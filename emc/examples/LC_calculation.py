from emc.mom2d import *
t=18e-6
w=890e-6
h1=144.6e-6
h2=200e-6
h3=144.6e-6
h_mask=30.e-6
er1=4.5
er2=4.7
er3=4.5
er_mask=3.5
s1=0.5e-3
s2=2*s1+w
d1 = w*5.0
all_W=2*d1+2*s1+w
d2=all_W-2.0*s2

isgrounded=False

config=Board()
config.layer(h1,er1)
config.conductor(1.e-6,d2,t,grounded=isgrounded)
config.conductor(s2,d2,t,grounded=isgrounded)

config.layer(h2,er2)
config.conductor(1e-6,d2,t,grounded=isgrounded)
config.conductor(s2,d2,t,grounded=isgrounded)

config.layer(h3,er3)
config.conductor(1e-6,d1,t,grounded=isgrounded)
config.conductor(s1,w,t)
config.conductor(s1,d1,t,grounded=isgrounded)

#config.conductor(1.e-6+d1+s1,w,t,grounded=isgrounded)

config.cover(h_mask,er_mask)


config.board2conf()
config.structure.set_subintervals(400)
print len(config.structure.list_diel)+len(config.structure.list_cond)

lc=RLGC(config.structure)
lc.calcLC()
if not isgrounded:
    m,n=lc.mL.shape
    print lc.mC[m-2,n-2]
    print lc.mL[m-2,n-2]
else:
   print lc.mC
   print lc.mL

