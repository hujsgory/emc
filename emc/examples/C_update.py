from emc.mom2d import *
t=18e-6
w=890e-6
h=144.6e-6
er=4.5
s=0.5e-3
d = w*5.0
all_W=2*d+2*s+w

isgrounded=False

board=Board()
board.layer(h,er)
board.conductor(d,w,t)
structure=board.board2structure()
structure.set_subintervals(100)

lc=RLGC(structure)
lc.calcC()
print lc.mC

er=5.5
board=Board()
board.layer(h,er)
board.conductor(d,w,t)
structure=board.board2structure()
structure.set_subintervals(100)

lc.update(structure)
print lc.mC