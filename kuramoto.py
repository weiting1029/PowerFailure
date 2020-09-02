#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 01:41:37 2020

@author: weiting
"""

# %% loading packages 
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import math 
import matplotlib.pyplot as plt
from multiprocessing import Pool 

##%%
## parameters setting one
#K = 6
#M = np.array([0.25, 0.75, 1, 1.25, 1, 0.75, 0.25])
#D = np.array([1, 2, 3, 2, 1, 2, 3])
#Ome = np.array([-3, -2, -1, 0.5, 1.5, 2.5, 3.5])
#n = np.size(M)
#pi = math.pi
#A = np.array([[0,0,0,0.95,0.5,0.2,0],[0,0,0,0,0.1,0,0.55],[0,0,0,0.55,0,0.3,0],[0.95,0,0.55,0,0,0.15,0],
#              [0.5,0.1,0,0,0,0,0],[0.2,0,0.3,0.15,0,0,0.9],[0,0.55,0,0,0,0.9,0]])
#t = 15
#nn = 100
#
#
## initial condtions
#theta0 = np.array([1, 1, -1, 1, 1, 1, -1])
#omega0 = np.array([-pi/5, pi/2, pi/5, 2*pi/3, -2*pi/5, -pi/3, -pi/7])
#
#
##replace the diagonal elements of A with zeros
#exA = A
#np.fill_diagonal(exA,0)
#
##calculate the sychronization frequency
#Omee =  np.sum(Ome)/np.sum(D)
#Ome0 = Ome-D*Omee
#




# %% parameter setting two


K = 1
n = 10
M = np.array([0.2228, 0.1607, 0.1899, 0.1517, 0.1379, 0.1846, 0.1401, 0.1289, 0.183, 2.6526])
D = np.array([0.0332,0.076,0.0862,0.0838,0.0674,0.0862,0.0743,0.0716,0.1101,0.1333])
Ome = np.zeros(10)


pi = math.pi


A = np.array([[0,2.484, 2.858, 2.372, 1.031, 0.222,2.13,13.057,3.676,7.834],
              [2.484,0,10.494,1.781,0.774,0.166,1.599,1.473,0.903,5.885],
              [2.858,10.494,0,2.522,1.096,0.236,2.264,1.729,1.144,4.772],
              [2.372,1.781,2.522,0,17.136,0.535,5.14,1.699,1.755,1.367],
              [1.031,0.774, 1.096,17.136,0 ,0.232, 2.235, 0.739, 0.763, 0.594],
              [0.222,0.166,0.236,0.535,0.232,0,1.574, 0.159, 0.164, 0.128],
              [2.13,1.599,2.264, 5.14, 2.235, 1.574, 0, 1.526, 1.576, 1.227],
              [13.057,1.473,1.729, 1.699,0.739, 0.159, 1.526, 0, 4.283, 4.159],
              [3.676, 0.903, 1.144, 1.755, 0.763, 0.164, 1.576, 4.283, 0, 1.37],
              [7.834, 5.885, 4.772, 1.367, 0.594, 0.128, 1.227, 4.159, 1.37 ,0]])

 
 
t = 15
nn = 100

theta0 = np.zeros(n)
omega0 = np.array([0.666, 1.484, -0.051, -0.153, -0.357, -0.201, 0.255, 0.465, -0.941, -0.012])




#replace the diagonal elements of A with zeros
exA = A
np.fill_diagonal(exA,0)

#calculate the sychronization frequency
Omee =  np.sum(Ome)/np.sum(D)
Ome0 = Ome-D*Omee









# %% test1
# def dU_dx(U, x):
#     # Here U is a vector such that y=U[0] and z=U[1]. This function should return [y', z']
#     return [U[1], -2*U[1] - 2*U[0] + np.cos(2*x)]
# U0 = [0, 0]
# xs = np.linspace(0, 10, 200)
# Us = odeint(dU_dx, U0, xs)
# ys = Us[:,0]

# plt.xlabel("x")
# plt.ylabel("y")
# plt.title("Damped harmonic oscillator")
# plt.plot(xs,ys);

# a,b,c,d = 1,1,1,1

# def dP_dt(P, t):
#     return [P[0]*(a - b*P[1]), -P[1]*(c - d*P[0])]

# ts = np.linspace(0, 12, 100)
# P0 = [1.5, 1.0]
# Ps = odeint(dP_dt, P0, ts)
# prey = Ps[:,0]
# predators = Ps[:,1]
# def kuramoto2nd(X, t):
#     dxdt = X[1]
#     dydt = 1 / m * (-d * X[1] + w - K * np.sum())





#%% define the solver 
def single(theta_j,omega_j,m,d,w,a_exj,theta_exj):
    dtheta_j = omega_j
    vec_temp = np.ones(n)*theta_j
    temp0 = np.sin(vec_temp-theta_exj)
    temp1 = np.multiply(a_exj,temp0)
    domega_j = 1/m*(-d*omega_j+w-K*np.sum(temp1))
    return[dtheta_j, domega_j]










def kuramoto2nd(X,t):
    theta = X[0:n]
    theta = theta-Omee*t
    omega = X[n:2*n]
    Y = np.zeros((n,2))
    for j in range(n):
        Y[j,:] = single(theta[j],omega[j],M[j],D[j],Ome0[j],exA[j,:],theta)
    
    return np.reshape(np.transpose(Y),2*n)




def kuramoto2nd2(X,t):
    theta = X[0:n]
    theta = theta-Omee*t
    omega = X[n:2*n]
    dtheta = omega 
    matrix1 = np.repeat(np.reshape(theta,(1,n)),n,axis=0)
    matrix2 = np.transpose(matrix1)-matrix1
    sinmatrix  = np.sin(matrix2)
    domega = (1/M)*(-D*omega+Ome0-K*np.sum(np.multiply(sinmatrix, exA),axis=1))
    return np.append(dtheta,domega)
    
                 
                  
                  
    
    
    
    
    
    
    
    
    




def solkuramoto(sol0,dt):
    return odeint(kuramoto2nd,sol0,dt)


def solkuramoto2(sol0,dt):
    return odeint(kuramoto2nd2,sol0,dt)

    



#%% a paralleized version 













# %% calculate the 1st derivative of omega
# domega = np.zeros((nn,n))
# for i in range(nn):
#     theta = sol_theta[i]
#     theta = theta-Omee*(i+1)*(t/nn)
#     omega = sol_omega[i]
#     Y = np.zeros((n,2))
#     for j in range(n):
#         Y[j,:] = single(theta[j],omega[j],M[j],D[j],Ome0[j],exA[j,:],theta)
    
#     domega[i] = np.transpose(Y[:,1])
    













 






















