'''
Created on 30 dic. 2018

@author: gsenno
'''

import qutip as qt
import numpy as np
import unittest
from main import constructJMatrix

class TestBellPolytope(unittest.TestCase):

    def testConstructJMatrix3inputs3outputs(self):
        N=3
        K=3
        d=2
        functional=np.random.uniform(-10,10,size=(N*K)**2)
        
        rho = np.random.uniform(-10,10,size=((d**2,d**2)))
        J=constructJMatrix(functional, N, K, d, rho)
        
        self.assertTrue(np.array_equal(J[0][0],-functional[0]*rho),"J[0][0] should be: - Bell_coeff_for_p(0,0|0,0)*rho")
        self.assertTrue(np.array_equal(J[1][2],-functional[5]*rho),"J[1][2] should be: - Bell_coeff_for_p(1,2|0,0)*rho")
        self.assertTrue(np.array_equal(J[0][3],-functional[9]*rho),"J[0][3] should be: - Bell_coeff_for_p(0,0|0,1)*rho")
    
    def testConstructJMatrix4inputs3outputs(self):
        N=4
        K=3
        d=2
        functional=np.random.uniform(-10,10,size=(N*K)**2)
        
        rho = np.random.uniform(-10,10,size=((d**2,d**2)))
        J=constructJMatrix(functional, N, K, d, rho)
        
        self.assertTrue(np.array_equal(J[0][0],-functional[0]*rho),"J[0][0] should be: - Bell_coeff_for_p(0,0|0,0)*rho")
        self.assertTrue(np.array_equal(J[1][2],-functional[5]*rho),"J[1][2] should be: - Bell_coeff_for_p(1,2|0,0)*rho")
        self.assertTrue(np.array_equal(J[0][3],-functional[9]*rho),"J[0][3] should be: - Bell_coeff_for_p(0,0|0,1)*rho")
    
    def testConstructJMatrix4inputs4outputs(self):
        N=4
        K=4
        d=2
        functional=np.random.uniform(-10,10,size=(N*K)**2)
        
        rho = np.random.uniform(-10,10,size=((d**2,d**2)))
        J=constructJMatrix(functional, N, K, d, rho)
    
        self.assertTrue(np.array_equal(J[0][0],-functional[0]*rho),"J[0][0] should be: - Bell_coeff_for_p(0,0|0,0)*rho")
        self.assertTrue(np.array_equal(J[1][2],-functional[6]*rho),"J[1][2] should be: - Bell_coeff_for_p(1,2|0,0)*rho")
        self.assertTrue(np.array_equal(J[0][4],-functional[16]*rho),"J[0][3] should be: - Bell_coeff_for_p(0,0|0,1)*rho")
# de

# def printMatrices(prob):
#     print('X00 = ', prob.matsolX(0))
#     print('X01 = ', prob.matsolX(1))
#     print('X02 = ', prob.matsolX(2))
#     print('X00+X01+X02 = ',prob.matsolX(0)+prob.matsolX(1)+prob.matsolX(2))
# 
#     print('-------------------------------------------------------------')    
#     
#     print('X10 = ', prob.matsolX(3))
#     print('X11 = ', prob.matsolX(4))
#     print('X12 = ', prob.matsolX(5))
#     print('X10+X11+X12 = ',prob.matsolX(3)+prob.matsolX(4)+prob.matsolX(5))
# 
#     print('-------------------------------------------------------------')    
#     
#     print('X20 = ', prob.matsolX(6))
#     print('X21 = ', prob.matsolX(7))
#     print('X22 = ', prob.matsolX(8))
#     print('X20+X21+X22 = ',prob.matsolX(6)+prob.matsolX(7)+prob.matsolX(8))
#     
#     print('-------------------------------------------------------------')    
# 
#     print('Y00 = ', prob.matsolY(0))
#     print('Y01 = ', prob.matsolY(1))
#     print('Y02 = ', prob.matsolY(2))
#     print('Y00+Y01+Y02 = ',prob.matsolY(0)+prob.matsolY(1)+prob.matsolY(2))
# 
#     print('-------------------------------------------------------------')    
#     
#     print('Y10 = ', prob.matsolY(3))
#     print('Y11 = ', prob.matsolY(4))
#     print('Y12 = ', prob.matsolY(5))
#     print('Y10+Y11+Y12 = ',prob.matsolY(3)+prob.matsolY(4)+prob.matsolY(5))
# 
#     print('-------------------------------------------------------------')    
#     
#     print('Y20 = ', prob.matsolY(6))
#     print('Y21 = ', prob.matsolY(7))
#     print('Y22 = ', prob.matsolY(8))
#     print('Y20+Y21+Y22 = ',prob.matsolY(6)+prob.matsolY(7)+prob.matsolY(8))