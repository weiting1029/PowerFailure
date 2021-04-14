import seaborn as sns

from greedy_algorithm import greedy_algorithm
from powerNetwork import getUndGraph, kron_reduction
from powerNetwork import networkTransform, getBuses, getLines
from powerNetworkSolver import PowerNetworkSolver
import numexpr as ne

# sns.set_theme()
# importing the function for multi-edge removal
import numpy as np
from disturbancesGnr import normaldisturbances
import math
from pypower.api import case39
from numpy.random import seed

#
# HERE: reset number of vml-threads
ne.set_vml_num_threads(8)

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
t = 2
nn = 100
dt = np.linspace(0, t, nn + 1)

sigma = 0.01
theta0 = np.zeros(ngnr39)
omega0 = np.zeros(ngnr39)
OMEGA = np.zeros(ngnr39)
model39 = PowerNetworkSolver(theta0, omega0, redA39, ngnr39, D, M, K, OMEGA)

seed(100)
disturbances = normaldisturbances(ngnr39, 1, sigma)
sol0 = np.pad(disturbances[0], (ngnr39, 0), 'constant', constant_values=(0, 0))
single_sol = model39.solkuramoto(sol0, dt)
sol_domega = model39.getDotOmega(single_sol[:, :ngnr39], single_sol[:, ngnr39:], nn)

check_times = 100
KK = 100  # repetition times
thres = np.array([0.2, 2])  # thres1 is for omega, thres2 is for omega_dot

un_graph = unG39
int_theta = theta0
int_omega = omega0
max_itr = 5
type_rate = 1

new_graph = greedy_algorithm(un_graph, int_theta, int_omega, D, M, K, OMEGA, KK, check_times, sigma, thres, t, nn,
                             normaldisturbances, max_itr, type_rate)
