#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 21:39:39 2020

@author: weiting
"""



import numpy as np 
from numpy.random import multivariate_normal


#changed 
def normaldisturbances(n,k,sigma):
    """
    ...
    Parameters
    ----------
    n : int
        number of nodes 
    k : int 
        sample size
    sigma : float 
        a scalar to control the covariance matrix 

    Returns
    -------
    ...
        (k,n)

    """
    cov = np.eye(n)*sigma 
    mu = np.zeros(n)
    return multivariate_normal(mu,cov,k)
    

    
def violationcheck(x,checktimes,thres):
    n = np.shape(x)[1]
    tt = np.shape(x)[0]
    result = np.zeros(n)
    checkstep = tt//checktimes
    for i in range(n):
        tempx= np.abs(x[:,i])
        for j in range(checktimes):
            if(tempx[j*checkstep]>thres):
                result[i]=1
                continue
    return result





def globalcheck(xx,checktimes,thres,n,tt):
    """
    a function returns all check results 

    Parameters
    ----------
    xx : an array 
        numerical solutions of the ode
    checktimes : int
        number of checks 
    thres : a (2,) array
        thres = [thres1, thres2]
    n : int 
        number of nodes
    tt : int
        number of integration steps 
        
    Returns
    -------
    a 3*n array 

    """
    result = np.zeros((3,n))
    checkstep = tt//checktimes
    absxx = np.abs(xx)
    for i in range(n):
        for k in range(2):  
            for j in range(checktimes):
                if(absxx[1+j*checkstep,k*n+i]>thres[k]):
                    result[k,i]=1
                    break  
        result[2,i] = np.max(result[:,i])
     
    return result 







def globalchecksubset(xx,checktimes,thres,n,tt):
    """
    a function returns all check results 

    Parameters
    ----------
    xx : an array 
        numerical solutions of the ode
    checktimes : int
        number of checks 
    thres : a (2,) array
        thres = [thres1, thres2]
    n : int 
        number of nodes
    tt : int
        number of integration steps 
        
    Returns
    -------
    a 3*n array 

    """
    result = np.zeros((3,n))
    checkstep = tt//checktimes
    absxx = np.abs(xx)
    for i in range(n):
        for k in range(2):  
            subset = absxx[1::checkstep,k*n+i]
            if (np.max(subset)>=1):
                result[k,i]=1
            # for j in range(checktimes):
            #     if(absxx[1+j*checkstep,k*n+i]>thres[k]):
            #         result[k,i]=1
            #         break  
        result[2,i] = np.max(result[:,i])
     
    return result 



# def globalcheck(xx,checktimes,thres,n,tt):
#     """
#     a function returns all check results 

#     Parameters
#     ----------
#     xx : an array 
#         numerical solutions of the ode
#     checktimes : int
#         number of checks 
#     thres : a (2,) array
#         thres = [thres1, thres2]
#     n : int 
#         number of nodes
#     tt : int
#         number of integration steps 
        
#     Returns
#     -------
#     a 3*n array 

#     """
#     result = np.zeros((3,n))
#     checkstep = tt//checktimes
#     for i in range(n):
#         for k in range(2):  
#             for j in range(checktimes):
#                 if(np.abs(xx[1+j*checkstep,i])>thres[0]):
#                     result[0,i]=1
#                     if(np.abs(xx[1+j*checkstep,n+i])>thres[1]):
#                         break  
#         result[2,i] = np.max(result[:,i])
     
#     return result 





                
                
# def check_joint_rocof(check1,check2):
#     joint_check = np.outer(check1,check2)
#     return joint_check
    