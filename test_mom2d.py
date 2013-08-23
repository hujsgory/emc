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
    def test_SmnAny2D(self):
        self.smn.SmnAny2D()
        self.assertEqual((numpy.around(self.smn.matrix,decimals=5)==[[3.84183,  1.05704],   [ 0.12258,  9.42478]]).all(),True)
    def test_SmnOrtho(self):
        self.smn.SmnOrtho()
        self.assertEqual((numpy.around(self.smn.matrix,decimals=5)==[[3.52549,  0.78204],   [-0.33587,  9.42478]]).all(),True)

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
        self.conf.diel(erp=2.0)
        self.conf.add(Section(Coord(0.0,0.0),Coord(0.0,1.0)))
        self.conf.add(Section(Coord(0.0,1.0),Coord(1.0,1.0)))
        self.conf.cond(erp=2.0)
        self.conf.add(Section(Coord(1.0,1.0),Coord(1.0,0.0)))
        self.conf.add(Section(Coord(1.0,0.0),Coord(0.5,0.5)))
        self.conf.cond(erp=3.0)
        self.conf.add(Section(Coord(0.5,0.5),Coord(0.0,0.0)))
        self.rlgc=RLCG(self.conf)
    def test_calcC(self):
        self.rlgc.calcC()
        self.assertEqual(self.rlgc.n_cond,2,self.rlgc.n_cond)
        self.assertEqual((numpy.around(self.rlgc.mC,5)==[[1.08443e-010,-3.39627e-011],[-4.60713e-011,1.25195e-010]]).all(),True)

unittest.main()
