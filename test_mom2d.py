#!/usr/bin/env python
# -*- coding: utf8 -*-
from mom2d import *
from math import *
import unittest
import numpy

class Test_mom2d_Coord(unittest.TestCase):
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

class Test_mom2d_Section(unittest.TestCase):
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

class Test_mom2d_Smn(unittest.TestCase):
    def setUp(self):
        self.s1=Section(Coord(0.0,0.0),Coord(1.0,1.0))
        self.s2=Section(Coord(1.0,1.0),Coord(2.0,0.0))
        self.conf=Conf()
        self.conf.add(self.s1, True , erp=2.0)
        self.conf.add(self.s2, False, erp=2.0)
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
    def test_Smn(self):
        self.assertEqual((numpy.around(SmnAny2D(self.conf),decimals=5)==[[3.84183,	1.05704],	[ 0.12258,	9.42478]]).all(),True)
        self.assertEqual((numpy.around(SmnOrtho(self.conf),decimals=5)==[[3.52549,	0.78204],	[-0.33587,	9.42478]]).all(),True)

class Test_mom2d_Conf(unittest.TestCase):
    def setUp(self):
        self.s1=Section(Coord(0.0,0.0),Coord(1.0,1.0))
        self.s2=Section(Coord(1.0,1.0),Coord(2.0,0.0))
        self.s3=Section(Coord(1.0,0.0),Coord(0.0,1.0))
        self.s4=Section(Coord(2.0,2.0),Coord(3.0,3.0))
        self.conf=Conf()
        self.conf.add(self.s1, True , erp=2.0)
        self.conf.add(self.s2, False, erp=2.0)
        self.conf.add(self.s4, True , erp=2.0)
    def test_intersection(self):
        self.assertRaises(ValueError,self.conf.add,self.s3,True,erp=2.0)

unittest.main()
