from numpy.random import multivariate_normal
import numpy as np
from numpy.random import seed
from numpy.linalg import matrix_power
from scipy.linalg import fractional_matrix_power


def normaldisturbances(n, k, sigma):
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
    seed(0)  # control the random disturbance generator
    cov = np.eye(n) * sigma
    mu = np.zeros(n)
    return multivariate_normal(mu, cov, k)


def correlated_disturbances(delta, alpha, tau, n, L, k):
    """
    :param delta: controls the global variance
    :param alpha: controls the covariance
    :param tau: smoothing variable
    :param n: the number of generators
    :param L: reduced Laplacian
    :param k: sample size
    :return: (k,n) random vectors
    """
    cov_matrix = delta * tau ** (2 * alpha) * fractional_matrix_power(L + (tau ** 2) * np.eye(n), -alpha)
    seed(100)
    mu = np.zeros(n)
    return multivariate_normal(mu, cov_matrix, k)


# import numpy as np

def correlation_from_covariance(covariance):
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    correlation[covariance == 0] = 0
    return correlation
