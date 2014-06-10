from emc.mom2d import *
t=18e-6
w=890e-6
h=144.6e-6
er=3.5
s=0.5e-3
d = w*5.0
all_W=2*d+2*s+w

board = Board()
board.layer(h, er)
board.conductor(d, w, t)
board.conductor(s, w, t)
struct = board.postprocess()
struct.set_subintervals(400)
lc=RLGC(struct)
lc.C()

w=1.5e-3
er=5.5
board = Board()
board.layer(h, er)
board.conductor(d, w, t)
board.conductor(s, w, t)
struct = board.postprocess()
struct.set_subintervals(400)
lc.update(struct)
print "iterative: \n", lc.mC
lc.C()
print "direct: \n", lc.mC
