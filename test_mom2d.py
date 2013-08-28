#!/usr/bin/env python
# -*- coding: utf8 -*-
from mom2d import *
from math import *
import unittest
import numpy

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
        self.assertTrue(((self.smn.matrix_S-[[3.84183,1.05704],[0.12258,9.42478]])<1e-5).all())
    def test_SmnOrtho(self):
        self.smn.SmnOrtho()
        self.assertTrue(((self.smn.matrix_S-[[3.52549,0.78204],[-0.33587,9.42478]])<1e-5).all())

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
        self.rlgc=RLCG(self.conf)
        self.rlgc.calcC()
    @unittest.skip('Variable n_cond is local in calc_C')
    def test_calcC1(self):
        self.assertEqual(self.rlgc.n_cond,2,self.rlgc.n_cond)
    @unittest.skip('Test for initial matrix_Q filling')
    def test_calcC2(self):
        self.assertTrue((self.rlgc.matrix_Q==[[Coef_C,0],[Coef_C,0],[Coef_C,0],[Coef_C,0],[Coef_C,0],[0,Coef_C],[0,Coef_C],[0,Coef_C],[0,Coef_C],[0,0],[0,0],[0,0]]).all())
    def test_calcC3(self):
        err=abs(self.rlgc.matrix_S-map(lambda x: map(float,x.split()),open('test_mom2d_test_RLGC_test_CalcC3.txt').readlines()))
        self.assertTrue((err<2e-4).all())
    def test_calcC4(self):
        self.assertTrue(((self.rlgc.mC-[[1.53885e-010,-4.32953e-011],[-6.27153e-011,2.06411e-010]])<1e-15).all())
    def test_calcL1(self):
        self.assertTrue(((self.rlgc.mL-[[6.35829e-007, 1.22435e-007],[ 1.28342e-007,7.29938e-007]])<1e-12).all())
unittest.main()
