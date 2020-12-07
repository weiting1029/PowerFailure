
import numpy as np
import copy
from numpy.linalg import inv 
import networkx as nx
from numpy.random import multivariate_normalx
from Disturbances import normaldisturbances


class PowerNetwork(object):
    def __init__(self, int_theta, int_omega, G, ngnr, D, M):
        super(PowerNetwork, self).__init__()
        self.int_omega = int_omega
        self.int_theta = int_theta
        self.G = G 
        self.ngnr= ngnr
        self.D = D 
        self.M = M 
        
    def get_num_nodes(self):
        return self.G.number_of_nodes()

    def kuramoto2nd(self,X,t):
    	Ome0 =  self.int_omega-self.D*(np.sum(self.int_omega)/np.sum(self.D))
    	theta = X[0:n]
    	theta = theta-np.sum(self.int_omega)/np.sum(self.D)*t
    	omega = X[n:2*n]
    	dtheta = omega 
    	matrix1 = np.repeat(np.reshape(theta,(1,n)),n,axis=0)
    	matrix1trix2 = np.transpose(matrix1)-matrix1
    	sinmatrix  = np.sin(matrix2)
    	domega = (1/M)*(-D*omega+Ome0-K*np.sum(np.multiply(sinmatrix, A),axis=1))
    	return np.append(dtheta,domega)

    def solkuramoto(self, sol0,dt):
    	return odeint(self.kuramoto2nd(),sol0,dt)



    def edge_removing(self, re_edge):
        node1 = re_edge[0]
        node2 = re_edge[1]
        copyunG = copy.deepcopy(self.G)
        copyunG.remove_edge(node1,node2)
        if(len(list(nx.connected_components(copyunG)))==1):
         	self.G.remove_edge(node1,node2)
        else:
        	print("The graph is not connected without edge ("+ str(node1) + ", " + str(node2) + ")")

    def kron_reduction(self):
     	n = self.get_num_nodes()
     	A = nx.to_numpy_matrix(self.G)
     	L = np.diag(np.squeeze(np.asarray(np.sum(A,axis = 0))))-A
     	redL = L[(n-self.ngnr):n,(n-self.ngnr):n] - np.transpose(L[:n-self.ngnr,n-self.ngnr:n])@inv(L[:n-self.ngnr,:n-self.ngnr])@L[:n-self.ngnr,n-self.ngnr:n]
     	redA = -redL
     	np.fill_diagonal(redA, 0)
     	return A, redL, redA

    def Simulation(self,KK,check_times,sigma,thres,t,nn,dtb_gnr):
    	dt = np.linspace(0, t, nn+1)
    	vcheck_omega = np.zeros(n)
		vcheck_theta = np.zeros(n)
		vcheck_any = np.zeros(n)
		mcheck_omega = np.zeros((n,n))
		mcheck_theta = np.zeros((n,n))
		mcheck_any = np.zeros((n,n))

		disturbances = dtb_gnr(n,KK,sigma)

		for k in range(KK):
    		sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
    		vec_sol = solkuramoto(sol0,dt)
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



		return {'vcheck_omega': vcheck_omega, 'vcheck_theta': vcheck_theta, 'vcheck_any': vcheck_any, 'mcheck_theta':mcheck_theta, 'mcheck_any':mcheck_any}







    





     
     












    


# 		