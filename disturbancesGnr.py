from numpy.random import multivariate_normal
import numpy as np
from numpy.random import seed
from numpy.linalg import matrix_power
from scipy.linalg import fractional_matrix_power

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
    seed(100)#control the random disturbance generator
    cov = np.eye(n)*sigma 
    mu = np.zeros(n)
    return multivariate_normal(mu,cov,k)



def correlated_disturbances(delta,alpha,tau,n,L,k):
    cov_matrix = delta*tau**(2*alpha)*fractional_matrix_power(L+tau**(2*alpha)*np.eye(n),-alpha)
    seed(100)
    mu = np.zeros
    return multivariate_normal(mu,cov_matrix,k)
