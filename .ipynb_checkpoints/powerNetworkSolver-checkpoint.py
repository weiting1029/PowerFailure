
import numpy as np
import copy
from numpy.linalg import inv 
import networkx as nx
from numpy.random import multivariate_normal
from disturbancesGnr import normaldisturbances
from scipy.integrate import odeint
from violationChecking import violationcheck,globalcheck,globalchecksubset


class PowerNetworkSolver(object):
	def __init__(self, int_theta, int_omega, A, ngnr, D, M, K):
		super(PowerNetworkSolver, self).__init__()
		self.int_omega = int_omega
		self.int_theta = int_theta
		self.A = A 
		self.ngnr= ngnr
		self.D = D 
		self.M = M
		self.K = K


	def kuramoto2nd(self,X,t):
		n = self.ngnr
		Ome0 =  self.int_omega-self.D*(np.sum(self.int_omega)/np.sum(self.D))
		theta = X[0:n]
		theta = theta-np.sum(self.int_omega)/np.sum(self.D)*t
		omega = X[n:2*n]
		dtheta = omega 
		matrix1 = np.repeat(np.reshape(theta,(1,n)),n,axis=0)
		matrix2 = np.transpose(matrix1)-matrix1
		sinmatrix  = np.sin(matrix2)
		domega = (1/self.M)*(-self.D*omega+Ome0-self.K*np.sum(np.multiply(sinmatrix, self.A),axis=1))
		return np.append(dtheta,domega)

	def solkuramoto(self, sol0, dt):
		return odeint(self.kuramoto2nd,sol0,dt)





	def Simulation(self,KK,check_times,sigma,thres,t,nn,dtb_gnr):
        starttime = time.time()
		n = self.ngnr
		dt = np.linspace(0, t, nn+1)
		vcheck_omega = np.zeros(n)
		vcheck_theta = np.zeros(n)
		vcheck_any = np.zeros(n)
		mcheck_omega = np.zeros((n,n))
		mcheck_theta = np.zeros((n,n))
		mcheck_any = np.zeros((n,n))

		disturbances = dtb_gnr(n,KK,sigma)

		for k in range(KK):
			sol0 = np.pad(disturbances[k], (n,0), 'constant', constant_values=(0,0))
			vec_sol = self.solkuramoto(sol0,dt)
			all_check = globalcheck(vec_sol, check_times, thres, n, nn)
			vcheck_theta += all_check[0]
			vcheck_omega += all_check[1]
			vcheck_any += all_check[2]

			mcheck_theta += np.outer(all_check[0],all_check[0])
			mcheck_omega += np.outer(all_check[1],all_check[1])
			mcheck_any += np.outer(all_check[2],all_check[2])


		vcheck_omega = vcheck_omega/KK
		vcheck_theta = vcheck_theta/KK
		vcheck_any = vcheck_any/KK
		mcheck_theta = mcheck_theta/KK
		mcheck_omega = mcheck_omega/KK
        mcheck_any = mcheck_any /KK

        stoptime = time.time()
        totaltime = stoptime-starttime
        
		return {'vcheck_omega': vcheck_omega, 'vcheck_theta': vcheck_theta,'vcheck_any': vcheck_any,'mcheck_theta':   mcheck_theta,'mcheck_any': mcheck_any, 'totaltime': "{:2.2}sec".format(totaltime)}













	 
	 















	# 		