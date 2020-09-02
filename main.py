#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 22:03:25 2020

@author: weiting
"""


#%% loading packages 
import numpy as np
import scipy as sp
from scipy.integrate import odeint
import math 
import matplotlib.pyplot as plt
import math 


from kuramoto import solkuramoto, solkuramoto2
from Disturbances import normaldisturbances, violationcheck, globalcheck,globalchecksubset
import time

#%%parameter setting

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
omega0 = np.zeros(n)




#replace the diagonal elements of A with zeros
# exA = A
# np.fill_diagonal(exA,0)

#calculate the sychronization frequency
# Omee =  np.sum(Ome)/np.sum(D)
# Ome0 = Ome-D*Omee


# #%%no disturbances
# t=15
# nn=1000
# dt = np.linspace(0,t, nn+1)
# omega0 = np.array([0.666,1.481,-0.051,-0.153, -0.357, -0.201, 0.255, 0.465, -0.941, -0.012])

# sol0 = np.zeros((2*n))
# sol0[0:n] = theta0
# sol0[n:2*n]= omega0


# starttime = time.time()
# vec_sol = solkuramoto(sol0,dt)
# t1 = time.time()-starttime

# starttime = time.time()
# vec_sol2 = solkuramoto2(sol0, dt)
# t2 = time.time()-starttime


# sol_theta = vec_sol[:,0:n]
# sol_omega = vec_sol[:,n:2*n]




# sol_theta2 = vec_sol2[:,0:n]
# sol_omega2 = vec_sol2[:,n:2*n]




# #calculate the sychronization frequency
# Omee =  np.sum(Ome)/np.sum(D)
# Ome0 = Ome-D*Omee

# check_times = 100

# # check_omega = violationcheck(sol_omega,100,1)
# # check_theta = violationcheck(sol_theta,100,1)
# # check_thetaomega = np.outer(check_omega, check_theta)
# # np.where(check_thetaomega == 2, check_thetaomega, np.zeros(np.size(check_thetaomega)))
# # np.where(check_thetaomega == 0, check_thetaomega, np.ones(np.size(check_thetaomega)))


# all_check = globalcheck(vec_sol, check_times, np.array([1,1]), n, nn)




#%%with omega disturbances
starttime = time.time()
omega0 = np.zeros(n)
KK = 2#repetition times  
t = 5
nn = 1000
dt = np.linspace(0, t, nn+1)
sigma = 0.1
disturb_norm = normaldisturbances(n,KK,sigma)
Ome0 = np.zeros(n)



opt_checks = np.array([10,25,50,100,200,250,500])

num_checks = len(opt_checks)
vcheck_omega = np.zeros((num_checks,n))
vcheck_theta = np.zeros((num_checks,n))
vcheck_any = np.zeros((num_checks,n))


for k in range(KK):
    sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
    vec_sol = solkuramoto2(sol0,dt)
    for j in range(num_checks):
        check_times = opt_checks[j]
        all_check = globalcheck(vec_sol, check_times, np.array([0.35,0.5]), n, nn)
        vcheck_theta[j] += all_check[0]
        vcheck_omega[j] += all_check[1]
        vcheck_any[j] += all_check[2]

vcheck_omega = vcheck_omega/KK
vcheck_theta = vcheck_theta/KK
vcheck_any = vcheck_any/KK
    

stoptime = time.time()
totaltime = stoptime-starttime
print("Total time needed= {:2.2}sec".format(totaltime))



#%%with omega disturbances
starttime = time.time()
omega0 = np.zeros(n)
t = 5
nn = 1000
check_times = 100
dt = np.linspace(0, t, nn+1)


# check_omega = np.zeros(n)
# check_theta = np.zeros(n)
# check_thetaomega = np.zeros(n)


# sol0 = np.zeros((2*n))
# sol0[0:n] = theta0
# sol0[n:2*n]= omega0+disturb_norm[0]
# vec_sol = solkuramoto(sol0,dt)
# sol_theta = vec_sol[:,0:n]
# sol_omega = vec_sol[:,n:2*n]
# check_omega = violationcheck(sol_omega,100,1)+check_omega
# check_theta = violationcheck(sol_theta,100,0.25)+check_theta
# temp_check_thetaomega = check_omega+check_theta
# temp_check_thetaomega= np.where(temp_check_thetaomega == 2, temp_check_thetaomega, np.zeros(n))
# temp_check_thetaomega=np.where(temp_check_thetaomega == 0, temp_check_thetaomega, np.ones(n))
# check_thetaomega = check_thetaomega + temp_check_thetaomega

Ome0 = np.zeros(n)
sigma = 0.1
NN = 40
opt_size = np.linspace(100, 4000,NN)
vcheck_omega = np.zeros((NN,n))
vcheck_theta = np.zeros((NN,n))
vcheck_any = np.zeros((NN,n))







for j in range(NN):
    KK = int(opt_size[j])
    disturb_norm = normaldisturbances(n,KK,sigma)
    # check_omega = np.zeros(n)
    # check_theta = np.zeros(n)
    # check_jointomega = np.zeros((n,n))
    # check_jointtheta = np.zeros((n,n))
    
    for k in range(KK):
        sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
        vec_sol = solkuramoto2(sol0,dt)
        all_check = globalcheck(vec_sol, check_times, np.array([0.35,0.5]), n, nn)
        vcheck_theta[j] += all_check[0]
        vcheck_omega[j] += all_check[1]
        vcheck_any[j] += all_check[2]
   
    vcheck_omega[j] = vcheck_omega[j]/KK
    vcheck_theta[j] = vcheck_theta[j]/KK
    vcheck_any[j]= vcheck_any[j]/KK
    
        
stoptime = time.time()
totaltime = stoptime-starttime

    
    
#%% sample size  = 1000, check times =  100
starttime = time.time()
omega0 = np.zeros(n)
Ome0 = np.zeros(n)
KK =1000#repetition times  
t = 5
nn = 1000
dt = np.linspace(0, t, nn+1)
sigma = 0.5
check_times = 100
thres = np.array([0.35,0.5])


vcheck_omega = np.zeros(n)
vcheck_theta = np.zeros(n)
vcheck_any = np.zeros(n)
mcheck_omega = np.zeros((n,n))
mcheck_theta = np.zeros((n,n))
mcheck_any = np.zeros((n,n))   


disturb_norm = normaldisturbances(n,KK,sigma)

for k in range(KK):
    sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
    vec_sol = solkuramoto2(sol0,dt)
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
print("Total time needed= {:2.2}sec".format(totaltime))





# #%%%test the subset method 

# starttime = time.time()
# omega0 = np.zeros(n)
# KK =1000#repetition times  
# t = 5
# nn = 1000
# dt = np.linspace(0, t, nn+1)
# sigma = 0.1
# disturb_norm = normaldisturbances(n,KK,sigma)
# Ome0 = np.zeros(n)
# check_times = 100
# vcheck_omega = np.zeros(n)
# vcheck_theta = np.zeros(n)
# vcheck_any = np.zeros(n)
# mcheck_omega = np.zeros((n,n))
# mcheck_theta = np.zeros((n,n))
# mcheck_any = np.zeros((n,n))   

# for k in range(KK):
#     sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
#     vec_sol = solkuramoto2(sol0,dt)
#     all_check = globalchecksubset(vec_sol, check_times, np.array([0.35,0.5]), n, nn)
#     vcheck_theta += all_check[0]
#     vcheck_omega += all_check[1]
#     vcheck_any += all_check[2]
    
#     mcheck_theta += np.outer(all_check[0],all_check[0])
#     mcheck_omega += np.outer(all_check[1],all_check[1])
#     mcheck_any += np.outer(all_check[2],all_check[2])
    
# vcheck_omega = vcheck_omega/KK
# vcheck_theta = vcheck_theta/KK
# vcheck_any = vcheck_any/KK
# mcheck_theta = mcheck_theta/KK
# mcheck_omega = mcheck_omega/KK
# mcheck_any = mcheck_any /KK


# stoptime = time.time()
# totaltime = stoptime-starttime
# print("Total time needed= {:2.2}sec".format(totaltime))






#%%%against different sigma  or thers
starttime = time.time()
omega0 = np.zeros(n)
KK =1000#repetition times  
t = 5
nn = 1000
dt = np.linspace(0, t, nn+1)
Ome0 = np.zeros(n)
check_times = 100


def simulation(sigma,thres):
    result = np.zeros(4)
    disturb_norm = normaldisturbances(n,KK,sigma)
    for k in range(KK):
        sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
        vec_sol = solkuramoto2(sol0,dt)
        all_check = globalcheck(vec_sol, check_times, thres, n, nn)
        result[0] += np.max(all_check[0])
        result[1] += np.max(all_check[1])
        result[2] += np.max(all_check[2])
        result[3] += all_check[0,5]
        
    result[0] = result[0]/KK
    result[1] = result[1]/KK
    result[2] = result[2]/KK
    result[3] = result[3]/KK


    return result 
        

#%%%sigma
    
candidate_sigma = np.array([0.01,0.05,0.1,0.5,1])
thres = np.array([0.35,0.5])
rates_sigma = np.zeros((len(candidate_sigma),4))

for i in range(len(candidate_sigma)):
    rates_sigma[i] = simulation(candidate_sigma[i],thres)
    


# for i in range(4):
#     plt.plot(candidate_sigma,rates_sigma[:,i],'-',label = 'rate'+str(i+1))


plt.plot(candidate_sigma,rates_sigma[:,0],'-',label = 'RoCof rate of the system')
plt.plot(candidate_sigma,rates_sigma[:,1],'-',label = 'Absolute frequency violation rate of the system')
plt.plot(candidate_sigma,rates_sigma[:,2],'-',label = 'Any failure rate of the system')
plt.plot(candidate_sigma,rates_sigma[:,3],'-',label = 'RoCof rate of Node 6')

plt.xlabel("sigma")
plt.ylabel("Different failure rates")
plt.ylim(0, 1)
plt.title("Thres 1 = 0.35, Thres 2 = 0.5")
plt.legend();


#%%rocof threshold 1


candidate_thres1 = np.array([0.2,0.3,0.4,0.5,0.6])
rate_thres1 = np.zeros((len(candidate_sigma),4))
thres2 = 0.5
sigma = 0.05

disturb_norm = normaldisturbances(n,KK,sigma)
for k in range(KK):
       sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
       vec_sol = solkuramoto2(sol0,dt)
       for i in range(len(candidate_thres1)):
           all_check = globalcheck(vec_sol, check_times, np.array([candidate_thres1[i],thres2]), n, nn)
           rate_thres1[i,0] += np.max(all_check[0])
           rate_thres1[i,1] += np.max(all_check[1])
           rate_thres1[i,2] += np.max(all_check[2])
           rate_thres1[i,3] += all_check[0,5]

            

rate_thres1 =rate_thres1/KK
plt.plot(candidate_thres1,rate_thres1[:,0],'-',label = 'RoCof rate of the system')
plt.plot(candidate_thres1,rate_thres1[:,1],'-',label = 'Absolute frequency violation rate of the system')
plt.plot(candidate_thres1,rate_thres1[:,2],'-',label = 'Any failure rate of the system')
plt.plot(candidate_thres1,rate_thres1[:,3],'-',label = 'RoCof rate of Node 6')

plt.xlabel("Threshold 1")
plt.ylabel("Different failure rates")
plt.ylim(0, 1)
plt.title("Thres 2 = 0.5, sigma = 0.05")
plt.legend();



#%%% against thres 2



candidate_thres2 = np.array([0.4,0.5,0.6,0.7,0.8])
rate_thres2 = np.zeros((len(candidate_sigma),4))
thres1 = 0.5
sigma = 0.05

disturb_norm = normaldisturbances(n,KK,sigma)
for k in range(KK):
       sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
       vec_sol = solkuramoto2(sol0,dt)
       for i in range(len(candidate_thres2))                      :
           all_check = globalcheck(vec_sol, check_times, np.array([thres1,candidate_thres2[i]]), n, nn)
           rate_thres2[i,0] += np.max(all_check[0])
           rate_thres2[i,1] += np.max(all_check[1])
           rate_thres2[i,2] += np.max(all_check[2])
           rate_thres2[i,3] += all_check[0,5]
            

rate_thres2 =rate_thres2/KK
plt.plot(candidate_thres2,rate_thres2[:,0],'-',label = 'RoCof rate of the system')
plt.plot(candidate_thres2,rate_thres2[:,1],'-',label = 'Absolute frequency violation rate of the system')
plt.plot(candidate_thres2,rate_thres2[:,2],'-',label = 'Any failure rate of the system')
plt.plot(candidate_thres2,rate_thres2[:,3],'-',label = 'RoCof rate of Node 6')

plt.xlabel("Threshold 2")
plt.ylabel("Different failure rates")
plt.ylim(0, 1)
plt.title("Thres 1 = 0.5, sigma = 0.05")
plt.legend();


#%%with standard deviations 
starttime = time.time()
omega0 = np.zeros(n)
t = 5
nn = 1000
check_times = 100
dt = np.linspace(0, t, nn+1)



Ome0 = np.zeros(n)
sigma = 0.1
NN = 40
opt_size = np.linspace(100, 4000,NN)
vcheck_omega = np.zeros((NN,2))
vcheck_theta = np.zeros((NN,2))
vcheck_any = np.zeros((NN,2))





for i in range(NN):
    KK = int(opt_size[i])
    disturb_norm = normaldisturbances(n,KK,sigma)
    for k in range(KK):
        sol0 = np.pad(disturb_norm[k], (n,0), 'constant', constant_values=(0,0))
        vec_sol = solkuramoto2(sol0,dt)
        all_check = globalcheck(vec_sol, check_times, np.array([0.35,0.5]), n, nn)
        vcheck_theta[i,0] += np.max(all_check[0])
        vcheck_omega[i,0] += np.max(all_check[1])
        vcheck_any[i,0] += np.max(all_check[2])
  


   
   
    vcheck_theta[i,0] = vcheck_theta[i,0]/KK
    vcheck_theta[i,1] = vcheck_theta[i,0] - vcheck_theta[i,0]**2
    vcheck_omega[i,0] = vcheck_omega[i,0]/KK
    vcheck_omega[i,1] = vcheck_omega[i,0]- vcheck_omega[i,0]**2
    vcheck_any[i,0] = vcheck_any[i,0]/KK
    vcheck_any[i,1] = vcheck_any[i,0] - vcheck_any[i,0]**2






    
        
stoptime = time.time()
totaltime = stoptime-starttime






#%% plot against number of, any violations

c =1.96# at 95% confidence level


plt.plot(opt_size,vcheck_any[:,0],'-',label = 'the any violation rate of the system')


plt.fill_between(opt_size, 
                  vcheck_any[:,0]-c*np.sqrt(vcheck_any[:,1])/np.sqrt(opt_size),
                  vcheck_any[:,0]+c*np.sqrt(vcheck_any[:,1])/np.sqrt(opt_size),
                  color='green', alpha=.5)

plt.xlabel("Sample size")
plt.ylabel("Any-violation rate")
plt.ylim(0.5, 1)
plt.legend();



#%% plot against sample size, RoCof

c =1.96# at 95% confidence level


plt.plot(opt_size,vcheck_omega[:,0],'-',label = 'the RoCoF rate of the system')


plt.fill_between(opt_size, 
                  vcheck_omega[:,0]-c*np.sqrt(vcheck_omega[:,1])/np.sqrt(opt_size),
                  vcheck_omega[:,0]+c*np.sqrt(vcheck_omega[:,1])/np.sqrt(opt_size),
                  color='green', alpha=.5)

plt.xlabel("Sample size")
plt.ylabel("RoCoF rate")
plt.ylim(0.5, 1)
plt.legend();


#%% plot against sample size, RoCof

c =1.96# at 95% confidence level


plt.plot(opt_size,vcheck_theta[:,0],'-',label = 'the absolute frequency violation rate of the system')


plt.fill_between(opt_size, 
                  vcheck_theta[:,0]-c*np.sqrt(vcheck_theta[:,1])/np.sqrt(opt_size),
                  vcheck_theta[:,0]+c*np.sqrt(vcheck_theta[:,1])/np.sqrt(opt_size),
                  color='green', alpha=.5)

plt.xlabel("Sample size")
plt.ylabel("Absolute frequency violation rate")
plt.ylim(0.5, 1)
plt.legend();




# %% plotting Theta


for i in range(n):
    theta = vec_sol[:,i]
    plt.plot(dt, theta, "-", label="omega"+ str(i+1))


plt.xlabel("Time")
plt.ylabel("Omega")    
plt.ylim(-1, 1)
plt.legend();
plt.title('Natural rotation frequencies')

# %% plotting Omega


for i in range(n):
    omega = vec_sol[:,n+i]  
    plt.plot(dt, omega, "-", label="omega dot "+ str(i+1))


plt.xlabel("Time")

plt.ylabel("Omega dot")
plt.ylim(-1, 1)
plt.legend();
plt.title('Derivatives of natural rotation frequencies')



# # %% plot dOmega
# for i in range(n):
#     domega_i = domega[:,i]
#     plt.plot(dt[1:], domega_i, "-", label="domega"+ str(i+1))


# plt.xlabel("Time")
# plt.ylabel("dOmega")
# plt.legend();



#%% plot against number of checks 1



for i in range(n):
    plt.plot(opt_checks,vcheck_any[:,i],'-',label = 'node'+str(i+1))

plt.xlabel("Sample size")
plt.ylabel("Any-violation rate")
plt.ylim(0, 1)
plt.legend();


#%% plot against number of checks 2
for i in range(n):
    plt.plot(opt_checks,vcheck_omega[:,i],'-',label = 'node'+str(i+1))

plt.xlabel("Sample size")
plt.ylabel("RoCoF violation rate")
plt.ylim(0, 1)
plt.legend();


#%% plot against number of checks 3


for i in range(n):
    plt.plot(opt_checks,vcheck_any[:,i],'-',label = 'node'+str(i+1))

plt.xlabel("Sample size")
plt.ylabel("Absolute frequency violation rate")
plt.ylim(0, 1)
plt.legend();


#%%%














