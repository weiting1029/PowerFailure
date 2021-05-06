from greedy_algorithm import greedy_algorithm
from powerNetwork import getUndGraph, kron_reduction
from powerNetwork import networkTransform, getBuses, getLines
from powerNetworkSolver import PowerNetworkSolver
import numexpr as ne
from functools import partial
from itertools import repeat
# sns.set_theme()
# importing the function for multi-edge removal
import numpy as np
from disturbancesGnr import normaldisturbances
import math
from pypower.api import case39
from numpy.random import seed
import pandas as pd
import time
from multiprocessing import Pool

# HERE: reset number of vml-threads
ne.set_vml_num_threads(6)

network39, subnetwork39 = networkTransform(case39())
df_lines39 = getLines(subnetwork39)
df_buses39 = getBuses(subnetwork39)
n39, ngnr39, unG39 = getUndGraph(df_buses39, df_lines39, network39)
A39, redL39, redA39 = kron_reduction(n39, ngnr39, unG39)
node_list = np.arange(ngnr39) + 1

K = 1
M = np.array([0.2228, 0.1607, 0.1899, 0.1517, 0.1379, 0.1846, 0.1401, 0.1289, 0.183, 2.6526])
D = np.array([0.0332, 0.076, 0.0862, 0.0838, 0.0674, 0.0862, 0.0743, 0.0716, 0.1101, 0.1333])
# Ome = np.zeros(n)
pi = math.pi
# t = 2
# nn = 100
theta0 = np.zeros(ngnr39)
omega0 = np.zeros(ngnr39)
OMEGA = np.zeros(ngnr39)
model39 = PowerNetworkSolver(theta0, omega0, redA39, ngnr39, D, M, K, OMEGA)


# KK = 100  # repetition times


# test_rates39 = model39.Simulation(KK, check_times, sigma, thres, t, nn, disturbances)

def main():
    Ncpus = 5
    N = 600
    pool = Pool(Ncpus)
    check_times = 10
    thres = np.array([0.2, 2])
    t = 2
    nn = 100
    sigma = 0.01
    parallel_func = partial(model39.parallelized_Simulation, check_times, thres, t, nn, sigma)
    pool.map(parallel_func, [N]*Ncpus)
    pool.close()
    pool.join()


if __name__ == "__main__":
    starttime = time.time()
    main()
    stoptime = time.time()
    totaltime = stoptime - starttime
    print("{:2.2}sec".format(totaltime))
