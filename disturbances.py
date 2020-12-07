from numpy.random import multivariate_normalx

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