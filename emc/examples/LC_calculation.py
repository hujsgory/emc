from emc.mom2d import *
t = 18e-6
w = 890e-6
h1 = 144.6e-6
h2 = 200e-6
h3 = 144.6e-6
h_mask = 30.e-6
er1 = 4.5
er2 = 4.7
er3 = 4.5
er_mask = 3.5
s1 = 0.5e-3
s2 = 2*s1+w
d1 = w*5.0
all_W = 2*d1+2*s1+w
d2 = all_W-2.0*s2

isgrounded=False

# Make a new object Board
board=Board()

# Method layer set a new layer and its parameters: thickness,
# relative dielectric permittivity,
# tangent of loss (default 0) and 
# relative permiability (default 1.0)
board.layer(h1, er1)

# Method conductor set rectangular conductor parameters:
# space to edge (x=0) or a previous conductor, width,
# thickness, immersion depth of the conductor to an insulator,
# whether grounded conductor (boolean type, default False), 
# and how to calculate the distance to the conductor: to edge or
# to center (boolean type, default False)
board.conductor(1.e-6, d2, t, grounded=isgrounded)
board.conductor(s2, d2, t, grounded=isgrounded)

board.layer(h2, er2)
board.conductor(1e-6, d2, t, grounded=isgrounded)
board.conductor(s2, d2, t, grounded=isgrounded)

board.layer(h3, er3)
board.conductor(1e-6, d1, t, grounded=isgrounded)
board.conductor(s1, w, t)
board.conductor(s1, d1, t, grounded=isgrounded)

board.cover(h_mask, er_mask)

# Build structure
conf=board.to_structure()
conf.set_subintervals(100)
print len(conf.list_diel)+len(conf.list_cond)

lc=RLGC(conf)
lc.calc_LC()
print lc.mC[-2, -2]
print lc.mL[-2, -2]

