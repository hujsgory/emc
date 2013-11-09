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

board=Board()
board.layer(h1,er1)
board.conductor(1.e-6,d2,t,grounded=isgrounded)
board.conductor(s2,d2,t,grounded=isgrounded)

board.layer(h2,er2)
board.conductor(1e-6,d2,t,grounded=isgrounded)
board.conductor(s2,d2,t,grounded=isgrounded)

board.layer(h3,er3)
board.conductor(1e-6,d1,t,grounded=isgrounded)
board.conductor(s1,w,t)
board.conductor(s1,d1,t,grounded=isgrounded)

#config.conductor(1.e-6+d1+s1,w,t,grounded=isgrounded)

board.cover(h_mask,er_mask)


conf=board.board2structure()
conf.set_subintervals(100)
print len(conf.list_diel)+len(conf.list_cond)

lc=RLGC(conf)
lc.calcLC()
if not isgrounded:
    m,n=lc.mL.shape
    print lc.mC[m-2,n-2]
    print lc.mL[m-2,n-2]
else:
   print lc.mC
   print lc.mL

