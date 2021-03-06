#!/usr/bin/env python
# -*- coding: utf8 -*-
from emc.mom2d import *
from math import *
import unittest
import numpy
def read_matrix(fname):
    return map(lambda x: map(float,x.split()),open(fname).readlines())
def read_vector(fname):
    return map(float,open(fname).readline().split())

  
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
        self.conf=Structure()
        self.conf.diel(erp=2.0)
        self.conf.add(self.s1)
        self.conf.cond(erp=2.0)
        self.conf.add(self.s2)
        self.conf.set_subintervals(3)
        self.smn=Smn(self.conf)
        self.smn.fill_SC()
    def test_fillS00(self):
        self.assertTrue(numpy.allclose(self.smn.matrix_SC.A00, read_matrix('_mom2d_Smn_CalcS00.txt')))
    def test_fillS01(self):
        self.assertTrue(numpy.allclose(self.smn.matrix_SC.A01, read_matrix('_mom2d_Smn_CalcS01.txt')))
    def test_fillS10(self):
        self.assertTrue(numpy.allclose(self.smn.matrix_SC.A10, read_matrix('_mom2d_Smn_CalcS10.txt')))
    def test_fillS11(self):
        self.assertTrue(numpy.allclose(self.smn.matrix_SC.A11, read_matrix('_mom2d_Smn_CalcS11.txt')))


class Test_Conf(unittest.TestCase):
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
        self.conf=Structure()
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


class Test_RLGC(unittest.TestCase):
    def setUp(self):
        self.conf=Structure()
        self.conf.diel(erp=2.0,mup=5.0)
        self.conf.add(Section(Coord(0.0,0.0),Coord(0.0,1.0)),1)
        self.conf.add(Section(Coord(0.0,1.0),Coord(1.0,1.0)),2)
        self.conf.cond(erp=2.0,mup=5.0)
        self.conf.add(Section(Coord(1.0,1.0),Coord(1.0,0.0)),2)
        self.conf.add(Section(Coord(1.0,0.0),Coord(0.5,0.5)),3)
        self.conf.cond(erp=3.0,mup=5.0)
        self.conf.add(Section(Coord(0.5,0.5),Coord(0.0,0.0)),4)
        self.rlgc=RLGC(self.conf)
        self.rlgc.calc_LC()
    def test_calcC1(self):
        self.assertEqual(self.rlgc.smn.n_cond,2,self.rlgc.smn.n_cond)
    def test_calcC2(self):
        self.assertTrue(numpy.allclose(self.rlgc.matrix_QC[0:9].T, read_matrix('_mom2d_RLGC_CalcC2.txt'), rtol=1e-6))
    def test_calcC5(self):
        self.assertTrue(numpy.allclose(self.rlgc.mC, read_matrix('_mom2d_RLGC_CalcC5.txt'), rtol=1e-6))
    def test_calcL1(self):
        self.assertTrue(numpy.allclose(self.rlgc.mL, read_matrix('_mom2d_RLGC_CalcL1.txt'), rtol=1e-6))
    def test_calcL3(self):
        self.assertTrue(numpy.allclose(self.rlgc.matrix_QL[0:9].T, read_matrix('_mom2d_RLGC_CalcL3.txt'), rtol=1e-6))


class Test_Board1(unittest.TestCase):
    def setUp(self):
        self.board=Board()
        self.board.layers=[{'height':1e-3,'er':2.0,'td':0.0,'mu':1.0,'cover':False,'cond':[{'space':1e-3,'width':0.25e-3,'thickness':100e-6,'depth':50e-6}]}]
    def test_layer1(self):
        self.assertRaises(ValueError,self.board.layer,height=50e-6,er=3.0)
    def test_conductor1(self):
        self.board.conductor(3e-3,1e-3,400e-6)
        self.assertEqual(self.board.layers[-1]['cond'][-1]['space'],3e-3)
        self.assertEqual(self.board.layers[-1]['cond'][-1]['width'],1e-3)
        self.assertEqual(self.board.layers[-1]['cond'][-1]['thickness'],400e-6)
        self.assertEqual(self.board.layers[-1]['cond'][-1]['depth'],0.0)
    def test_conductor2(self):
        self.assertRaises(ValueError,self.board.conductor,space=3e-3,width=1e-3,thickness=400e-6,depth=1e-3)
    def test_conductor3(self):
        self.assertRaises(ValueError,self.board.conductor,space=0.0,width=1e-3,thickness=400e-6,depth=0.0)
    def test_conductor4(self):
        self.assertRaises(ValueError,self.board.conductor,space=3e-3,width=0.0,thickness=400e-6,depth=0.0)
    def test_conductor5(self):
        self.assertRaises(ValueError,self.board.conductor,space=3e-3,width=1e-3,thickness=0.0,depth=0.0)
    def test_conductor6(self):
        self.assertRaises(ValueError,self.board.conductor,space=3e-3,width=1e-3,thickness=400e-6,depth=-0.1)
    @unittest.skip('Variable max_x is local')
    def test_board2conf(self):
        self.board.board2structure()
        self.assertAlmostEqual(self.board.max_x,2.25e-3,15)


class Test_Board2(unittest.TestCase):
    def setUp(self):
        self.board=Board()
        self.board.layer(990e-6,4.3)
        self.board.conductor(800e-6,401e-6,18e-6)
        self.board.medium['mu']=1.00000037
    def test_1(self):
        conf=self.board.to_structure()
        x0,x1,x2,x3=0.0,800e-6,1201e-6,2001e-6
        y1,y2=990e-6,1008e-6
        er1,er2=4.3,1.0
        td1,td2=0.0,0.0
        mu1,mu2=1.0,1.00000037
        answ1=[{'section':Section(Coord(x1,y1),Coord(x2,y1)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1},'obj_count':1,'sect_count':0,'grounded':False},\
               {'section':Section(Coord(x2,y1),Coord(x2,y2)),'n_subint':1,'mat_param':{'erp':er2,'tdp':td2,'mup':mu2},'obj_count':1,'sect_count':1,'grounded':False},\
               {'section':Section(Coord(x1,y2),Coord(x1,y1)),'n_subint':1,'mat_param':{'erp':er2,'tdp':td2,'mup':mu2},'obj_count':1,'sect_count':2,'grounded':False},\
               {'section':Section(Coord(x2,y2),Coord(x1,y2)),'n_subint':1,'mat_param':{'erp':er2,'tdp':td2,'mup':mu2},'obj_count':1,'sect_count':3,'grounded':False}]
        answ2=[{'section':Section(Coord(x0,y1),Coord(x1,y1)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':2,'sect_count':0},\
               {'section':Section(Coord(x2,y1),Coord(x3,y1)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':2,'sect_count':1}]
        self.assertTrue(conf.list_cond==answ1)
        self.assertTrue(conf.list_diel==answ2)
    def test_2(self):
        self.board.cover(31e-6,4.8)
        conf=self.board.to_structure()
        er1,er2=4.8,1.0
        td1,td2=0.0,0.0
        mu1,mu2=1.0,1.00000037
        x0,x1,x2,x3=0.0,769e-6,1232e-6,2001e-6
        y2,y3=1021e-6,1039e-6
        answ=[\
              {'section':Section(Coord(x0,y2),Coord(x1,y2)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':3,'sect_count':0},\
              {'section':Section(Coord(x1,y2),Coord(x1,y3)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':3,'sect_count':1},\
              {'section':Section(Coord(x1,y3),Coord(x2,y3)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':3,'sect_count':2},\
              {'section':Section(Coord(x2,y3),Coord(x2,y2)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':3,'sect_count':3},\
              {'section':Section(Coord(x2,y2),Coord(x3,y2)),'n_subint':1,'mat_param':{'erp':er1,'tdp':td1,'mup':mu1,'erm':er2,'tdm':td2,'mum':mu2},'obj_count':3,'sect_count':4},\
              ]
        self.assertTrue(conf.list_diel[2:]==answ)


class Test_Matrix(object):
    def SetUp(self):
        pass
    def test(self):
        pass

unittest.main()
