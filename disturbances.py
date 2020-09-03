#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug  2 21:39:39 2020

@author: weiting
"""



import numpy as np 
from numpy.random import multivariate_normal


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
    

 

