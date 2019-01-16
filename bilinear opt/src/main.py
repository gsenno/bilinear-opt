'''


Created on 23 Oct 2018





@author: gsenno


'''

import cvxopt as cvx
import numpy as np
import picos as sdp
from problem import *
from semidefinite import *
from bellpolytope import BellPolytope



def findMaxQuantumViolation(bellIneq,N,K,d,rho):
    prob = QuantumBilinearProblem()
    dims = [d*np.ones((1,N*K)),d*np.ones((1,N*K))]
    J = constructJMatrix(bellIneq,N,K,d,rho)
    prob.init_matrix_form(dims,J);
    
    
def constructJMatrix(bellIneq,N,K,d,rho):
    J = np.zeros((N*K,N*K,d*d,d*d))
    bellCoeff = 0;
    for x in range(N):
        for y in range(N):
            for a in range(K):
                for b in range(K):
                    J[a+K*x][b+K*y]=bellIneq[bellCoeff]*rho*(-1)
                    bellCoeff += 1
    return J.tolist()


def computeBellValue(N, K, rho, prob, functional):
    products = []
    for x in range(N):
        for y in range(N):
            for a in range(K):
                for b in range(K):
                    products.append(np.trace(rho*np.kron(prob.matsolX(a+N*x),prob.matsolY(b+N*y))))
    return np.dot(np.array(functional),products)


if __name__ == '__main__':
    parties = 2
    N = 4
    K = 3
    dim = 2
    
    poly=BellPolytope(N,K)   
    vertices=poly.getVertices()
    distributions=np.matrix(vertices)
    NumberOfCoefficients=(N*K)**2
             
    with open('results.txt','w') as f:
        f.write('Random functionals \n')
    
    rho=np.matrix([[0,0,0,0],[0,1/2,-1/2,0],[0,-1/2,1/2,0],[0,0,0,0]])
    #Ineqs = np.loadtxt( 'Ineqs.txt')
    for i in range (0,2): 
        functional=np.random.uniform(-N**2*(N**2-2),1,size=NumberOfCoefficients)
        values=np.dot(distributions,np.transpose(functional))
        c=np.amax(values)
        if c==0:
            functional=np.random.uniform(-N**2*(N**2-2),1,size=NumberOfCoefficients)
            values=np.dot(distributions,np.transpose(functional))
            c=np.amax(values)
        
        BellNormalisedFunctional=(functional/c).tolist()
        print(BellNormalisedFunctional[12])
    
        
        
        
        prob = QuantumBilinearProblem()
        #functional=Ineqs[i][:]
        #print(len(functional))
        J=constructJMatrix(BellNormalisedFunctional, N, K, dim, rho)
        dims = (dim*np.ones((parties,N*K),np.int8)).tolist()
        prob.init_matrix_form(dims,J)

        for j in range (0, (K*N)-1):
            prob.add_constraint(prob.matvarX(j) == prob.matvarX(j).H)
            prob.add_constraint(prob.matvarX(j) >> 0)
            prob.add_constraint(prob.matvarY(j) == prob.matvarY(j).H)
            prob.add_constraint(prob.matvarY(j) >> 0)
    
        '''prob.add_constraint(prob.matvarX(1) == prob.matvarX(1).H)
        prob.add_constraint(prob.matvarX(1) >> 0)
        prob.add_constraint(prob.matvarX(2) == prob.matvarX(2).H)
        prob.add_constraint(prob.matvarX(2) >> 0)
        prob.add_constraint(prob.matvarX(3) == prob.matvarX(3).H)
        prob.add_constraint(prob.matvarX(3) >> 0)
        prob.add_constraint(prob.matvarX(4) == prob.matvarX(4).H)
        prob.add_constraint(prob.matvarX(4) >> 0)
        prob.add_constraint(prob.matvarX(4)-10 ** (-3) >>0)
        prob.add_constraint(prob.matvarX(5) == prob.matvarX(5).H)
        prob.add_constraint(prob.matvarX(5) >> 0)
        prob.add_constraint(prob.matvarX(5)-10 ** (-3) >>0)
    
        prob.add_constraint(prob.matvarX(6) == prob.matvarX(6).H)
        prob.add_constraint(prob.matvarX(6) >> 0)
        prob.add_constraint(prob.matvarX(6)-10 ** (-3) >>0)
        prob.add_constraint(prob.matvarX(7) == prob.matvarX(7).H)
        prob.add_constraint(prob.matvarX(7) >> 0)
        prob.add_constraint(prob.matvarX(8) == prob.matvarX(8).H)
        prob.add_constraint(prob.matvarX(8) >> 0)
        prob.add_constraint(prob.matvarX(9) == prob.matvarX(9).H)
        prob.add_constraint(prob.matvarX(9) >> 0)
        prob.add_constraint(prob.matvarX(10) == prob.matvarX(10).H)
        prob.add_constraint(prob.matvarX(10) >> 0)
        prob.add_constraint(prob.matvarX(11) == prob.matvarX(11).H)
        prob.add_constraint(prob.matvarX(11) >> 0)'''

#         
        
        '''prob.add_constraint(prob.matvarY(1) == prob.matvarY(1).H)
        prob.add_constraint(prob.matvarY(1) >> 0)
        prob.add_constraint(prob.matvarY(2) == prob.matvarY(2).H)
        prob.add_constraint(prob.matvarY(2) >> 0)
        prob.add_constraint(prob.matvarY(2)-10 ** (-3) >>0)
        prob.add_constraint(prob.matvarY(3) == prob.matvarY(3).H)
        prob.add_constraint(prob.matvarY(3) >> 0)
        prob.add_constraint(prob.matvarY(4) == prob.matvarY(4).H)
        prob.add_constraint(prob.matvarY(4) >> 0)
        prob.add_constraint(prob.matvarY(5) == prob.matvarY(5).H)
        prob.add_constraint(prob.matvarY(5) >> 0)
        prob.add_constraint(prob.matvarY(6) == prob.matvarY(6).H)
        prob.add_constraint(prob.matvarY(6) >> 0)
        prob.add_constraint(prob.matvarY(7) == prob.matvarY(7).H)
        prob.add_constraint(prob.matvarY(7) >> 0)
        prob.add_constraint(prob.matvarY(8) == prob.matvarY(8).H)
        prob.add_constraint(prob.matvarY(8) >> 0)
        prob.add_constraint(prob.matvarY(8) == prob.matvarY(8).H)
        prob.add_constraint(prob.matvarY(8) >> 0)'''
    

        prob.add_constraint(prob.varX()[0]+prob.varX()[4] + prob.varX()[8] + prob.varX()[12]== 2./np.sqrt(2))
        prob.add_constraint(prob.varX()[1]+prob.varX()[5] + prob.varX()[9] + prob.varX()[13]< 10**-20.)
        prob.add_constraint(prob.varX()[1]+prob.varX()[5] + prob.varX()[9] + prob.varX()[13]> 0.)
        prob.add_constraint(prob.varX()[2]+prob.varX()[6] + prob.varX()[10] + prob.varX()[14]== 0.)
        prob.add_constraint(prob.varX()[3]+prob.varX()[7] + prob.varX()[11]  + prob.varX()[15]== 0.)

        prob.add_constraint(prob.varX()[16]+prob.varX()[20]+prob.varX()[24] + prob.varX()[28]== 2./np.sqrt(2))
        prob.add_constraint(prob.varX()[17]+prob.varX()[21]+prob.varX()[25] +prob.varX()[29]== 0)
        prob.add_constraint(prob.varX()[18]+prob.varX()[22]+prob.varX()[26] + prob.varX()[30]== 0.)
        prob.add_constraint(prob.varX()[19]+prob.varX()[23]+prob.varX()[27]  + prob.varX()[31]== 0.)

        prob.add_constraint(prob.varX()[32]+prob.varX()[36]+prob.varX()[40] + prob.varX()[44]== 2./np.sqrt(2))
        prob.add_constraint(prob.varX()[33]+prob.varX()[37]+prob.varX()[41] + prob.varX()[45]== 0)
        prob.add_constraint(prob.varX()[34]+prob.varX()[38]+prob.varX()[42] + prob.varX()[46]== 0.)
        prob.add_constraint(prob.varX()[35]+prob.varX()[39]+prob.varX()[43] + prob.varX()[47]== 0.)

        prob.add_constraint(prob.varY()[0]+prob.varY()[4] + prob.varY()[8] + prob.varY()[12]== 2./np.sqrt(2))
        prob.add_constraint(prob.varY()[1]+prob.varY()[5] + prob.varY()[9] + prob.varY()[13]< 10**-20.)
        prob.add_constraint(prob.varY()[1]+prob.varY()[5] + prob.varY()[9] + prob.varY()[13]> 0.)
        prob.add_constraint(prob.varY()[2]+prob.varY()[6] + prob.varY()[10] + prob.varY()[14]== 0.)
        prob.add_constraint(prob.varY()[3]+prob.varY()[7] + prob.varY()[11]  + prob.varY()[15]== 0.)

        prob.add_constraint(prob.varY()[16]+prob.varY()[20]+prob.varY()[24] + prob.varY()[28]== 2./np.sqrt(2))
        prob.add_constraint(prob.varY()[17]+prob.varY()[21]+prob.varY()[25] +prob.varY()[29]== 0)
        prob.add_constraint(prob.varY()[18]+prob.varY()[22]+prob.varY()[26] + prob.varY()[30]== 0.)
        prob.add_constraint(prob.varY()[19]+prob.varY()[23]+prob.varY()[27]  + prob.varY()[31]== 0.)

        prob.add_constraint(prob.varY()[32]+prob.varY()[36]+prob.varY()[40] + prob.varY()[44]== 2./np.sqrt(2))
        prob.add_constraint(prob.varY()[33]+prob.varY()[37]+prob.varY()[41] + prob.varY()[45]== 0)
        prob.add_constraint(prob.varY()[34]+prob.varY()[38]+prob.varY()[42] + prob.varY()[46]== 0.)
        prob.add_constraint(prob.varY()[35]+prob.varY()[39]+prob.varY()[43] + prob.varY()[47]== 0.)
    
    

        prob.solve(verb = 1, maxiter = STD_MAXITER, allowedgap = STD_ALLOWEDGAP)
    
        print("Bell value=",computeBellValue(N,K,rho,prob,BellNormalisedFunctional))
        if computeBellValue(N, K, rho, prob, BellNormalisedFunctional)<-2:
            with open('results.txt','a') as f:
                np.savetxt(f, BellNormalisedFunctional, fmt='%.2f')
        else:
            with open('results.txt','a') as f:
                f.write('Not this one \n')
#printMatrices(prob)
        for i in range(9):
            if np.all(np.linalg.eigvals(prob.matsolX(i)) >= 0)==False:
                print("X",i,prob.matsolX(i))
            if np.all(np.linalg.eigvals(prob.matsolY(i)) >= 0)==False:
                print("Y",i,prob.matsolY(i))
    