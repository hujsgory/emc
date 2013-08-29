#!/usr/bin/env python
# -*- coding: utf8 -*-
from mom2d import *
from math import *
import unittest
import numpy
def read_matrix(fname):
    return map(lambda x: map(float,x.split()),open('./test/'+fname).readlines())
def read_vector(fname):
    return map(float,open('./test/'+fname).readline().split())
        
class Test_Coord(unittest.TestCase):
    def setUp(self):
        self.m=Coord(1.0,1.0)
        self.n=Coord(0.0,0.0)
    def test_sint(self):
        self.assertEqual(self.m.sint(self.n)  ,-1.0/sqrt(2.0) , (self.m.sint(self.n),    -1.0/sqrt(2.0)))
        self.assertEqual(self.n.sint(self.m)  , 1.0/sqrt(2.0) , (self.n.sint(self.m),     1.0/sqrt(2.0)))
    def test_cost(self):
        self.assertEqual(self.m.cost(self.n)  ,-1.0/sqrt(2.0) , (self.m.cost(self.n),    -1.0/sqrt(2.0)))
        self.assertEqual(self.n.cost(self.m)  , 1.0/sqrt(2.0) , (self.n.cost(self.m),     1.0/sqrt(2.0)))
    def test_center(self):
        self.assertEqual(self.n.center(self.m), Coord(0.5,0.5), (self.n.center(self.m).x, self.n.center(self.m).y))
    def test_len(self):
        self.assertEqual(self.m.len(self.n)   , sqrt(2.0)     , (self.m.len(self.n),      sqrt(2.0)))

class Test_Section(unittest.TestCase):
    def setUp(self):
        self.s=Section(Coord(1.0,1.0),Coord(0.0,0.0))
    def test_center(self):
        self.assertEqual(self.s.center,   Coord(0.5,0.5), (self.s.center.x,self.s.center.y))
    def test_sint(self):
        self.assertEqual(self.s.sint,    -1.0/sqrt(2.0) , (self.s.sint,-1.0/sqrt(2.0)))
    def test_cost(self):
        self.assertEqual(self.s.cost,    -1.0/sqrt(2.0) , (self.s.cost,-1.0/sqrt(2.0)))
    def test_len(self):
        self.assertEqual(self.s.len,      sqrt(2.0)     , (self.s.len,sqrt(2.0)))
    def test_dx(self):
        self.assertEqual(self.s.dx,      -1.0           , self.s.dx)
    def test_dy(self):
        self.assertEqual(self.s.dy,      -1.0           , self.s.dy)
    def test_getSubinterval(self):
        self.s=Section(Coord(0.5,1.0),Coord(0.0,0.0))
        self.assertEqual(self.s.getSubinterval(2,4),Section(Coord(0.25,0.5),Coord(0.125,0.25)))
        self.assertRaises(ValueError,self.s.getSubinterval,4,4)

class Test_Smn(unittest.TestCase):
    def setUp(self):
        self.s1=Section(Coord(0.0,0.0),Coord(1.0,1.0))
        self.s2=Section(Coord(1.0,1.0),Coord(2.0,0.0))
        self.conf=Conf()
        self.conf.diel(erp=2.0)
        self.conf.add(self.s1)
        self.conf.cond(erp=2.0)
        self.conf.add(self.s2)
        self.smn=Smn(self.conf)
    def test_SmnAny2D(self):
        self.smn.SmnAny2D()
        self.assertTrue((abs(self.smn.matrix_S-[[3.84183,1.05704],[0.12258,0]])<1e-5).all())
        self.assertTrue(abs(self.smn.diag_S11_C[0]-9.42478)<1e-5)
    def test_SmnOrtho(self):
        self.smn.SmnOrtho()
        self.assertTrue((abs(self.smn.matrix_S-[[3.52549,0.78204],[-0.33587,9.42478]])<1e-5).all())

class Test_mom2d_Conf(unittest.TestCase):
    def setUp(self):
        self.c1=Coord(0.0,0.0)
        self.c5=Coord(1.0,1.0)
        self.c6=Coord(0.25,0.25)
        self.c7=Coord(0.5,0.5)
        self.cx2=Coord(0.5,0.0)
        self.cx3=Coord(1.0,0.0)
        self.cx4=Coord(2.0,0.0)
        self.cy2=Coord(0.0,0.5)
        self.cy3=Coord(0.0,1.0)
        self.cy4=Coord(0.0,2.0)
        self.conf=Conf()
        self.conf.cond(erm=2.0)
    def test_intersection_1(self):
        self.conf.add(Section(self.c1,self.cx3))
        self.conf.add(Section(self.cx3,self.cx4))
        self.conf.add(Section(self.cx3,self.c5))
        self.assertRaises(ValueError,self.conf.add,Section(self.cx2,self.cx3))
    def test_intersection_2(self):
        self.conf.add(Section(self.c1,self.cy3))
        self.conf.add(Section(self.cy3,self.cy4))
        self.conf.add(Section(self.cy3,self.c5))
        self.assertRaises(ValueError,self.conf.add,Section(self.cy2,self.cy3))
    def test_intersection_3(self):
        self.conf.add(Section(self.c6,self.c5))
        self.assertRaises(ValueError,self.conf.add,Section(self.c1,self.c7))

class Test_RLCG(unittest.TestCase):
    def setUp(self):
        self.conf=Conf()
        self.conf.diel(erp=2.0,mup=5.0)
        self.conf.add(Section(Coord(0.0,0.0),Coord(0.0,1.0)),1)
        self.conf.add(Section(Coord(0.0,1.0),Coord(1.0,1.0)),2)
        self.conf.cond(erp=2.0,mup=5.0)
        self.conf.add(Section(Coord(1.0,1.0),Coord(1.0,0.0)),2)
        self.conf.add(Section(Coord(1.0,0.0),Coord(0.5,0.5)),3)
        self.conf.cond(erp=3.0,mup=5.0)
        self.conf.add(Section(Coord(0.5,0.5),Coord(0.0,0.0)),4)
        self.rlcg=RLCG(self.conf)
        self.rlcg.calcLC()
    @unittest.skip('Variable n_cond is local in calc_C')
    def test_calcC1(self):
        self.assertEqual(self.rlgc.n_cond,2,self.rlgc.n_cond)
    def test_calcC2(self):
        err=abs(numpy.transpose(self.rlcg.matrix_QC[0:9,0:2])-read_matrix('_mom2d_RLCG_CalcC2.txt'))
        self.assertTrue((err<2e-15).all())
    def test_calcC3(self):
        err=abs(self.rlcg.matrix_S-read_matrix('_mom2d_RLCG_CalcC3.txt'))
        self.assertTrue((err<2e-4).all())
    def test_calcC4(self):
        err=abs(self.rlcg.diag_S11_C-read_vector('_mom2d_RLCG_CalcC4.txt'))
        self.assertTrue((err<2e-4).all())
    def test_calcC5(self):
        err=abs(self.rlcg.mC-read_matrix('_mom2d_RLCG_CalcC5.txt'))
        self.assertTrue((err<1e-15).all())
    def test_calcL1(self):
        err=abs(self.rlcg.mL-read_matrix('_mom2d_RLCG_CalcL1.txt'))
        self.assertTrue((err<1e-12).all())
    def test_calcL2(self):
        err=abs(self.rlcg.diag_S11_L-read_vector('_mom2d_RLCG_CalcL2.txt'))
        self.assertTrue((err<2e-4).all())
    def test_calcL3(self):
        err=abs(numpy.transpose(self.rlcg.matrix_QL[0:9,0:2])-read_matrix('_mom2d_RLCG_CalcL3.txt'))
        self.assertTrue((err<2e-16).all())
unittest.main()
