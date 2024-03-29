import pandas as pd
from powerNetworkSolver import PowerNetworkSolver
from powerNetwork import networkTransform, getGenerators, getBuses, getLines, getNumBuses, getNumLines
from powerNetwork import getUndGraph, kron_reduction, edge_removing
from disturbancesGnr import normaldisturbances
from violationChecking import globalcheck
import seaborn as sns
from disturbancesGnr import correlation_from_covariance

sns.set_theme()
# importing the function for multi-edge removal
from powerNetwork import multi_edge_removing
import matplotlib.pyplot as plt

import numpy as np
import copy
from numpy.linalg import inv
import networkx as nx
from numpy.random import multivariate_normal
from disturbancesGnr import normaldisturbances, correlated_disturbances
from scipy.integrate import odeint
from violationChecking import globalcheck
import time
import math
from functools import reduce

from pypower.api import case39
import numexpr as ne

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
t = 10
nn = 1000
dt = np.linspace(0, t, nn + 1)

sigma = 0.01
theta0 = np.zeros(ngnr39)
omega0 = np.zeros(ngnr39)
OMEGA = np.zeros(ngnr39)
model39 = PowerNetworkSolver(theta0, omega0, redA39, ngnr39, D, M, K, OMEGA)

from numpy.random import seed

# generate numerical tracks
# seed(100)
disturbances = correlated_disturbances(delta=1, alpha=1.5, tau=1.5, n=ngnr39, L=redL39, k=1)
# disturbances = normaldisturbances(n=ngnr39,k=1,sigma=sigma)
sol0 = np.pad(disturbances[0], (ngnr39, 0), 'constant', constant_values=(0, 0))
single_sol = model39.solkuramoto(sol0, dt)
sol_domega = model39.getDotOmega(single_sol[:, :ngnr39], single_sol[:, ngnr39:], nn)

# fig = plt.figure()
# for i in range(ngnr39):
#     theta = single_sol[:,i]
#     plt.plot(dt, theta, "-", label="Theta"+ str(i+1))
#
# plt.xlabel("Time")
# plt.ylabel("Theta")
# plt.ylim(-0.2, 1)
# plt.legend();
# plt.title('Phase angles')
# plt.show()
# plt.interactive(True)
# plt.savefig('theta_corcov.png')
#
#
fig = plt.figure()
for i in range(ngnr39):
    omega = single_sol[:,ngnr39+i]
    plt.plot(dt, omega, "-", label="Omega"+ str(i+1))


plt.xlabel("Time")

plt.ylabel("Omega")
plt.ylim(-2, 2)
plt.legend()
plt.title('Natural rotation frequencies')
# plt.savefig('omega_corcov.png')
#
#
#
#
# fig = plt.figure()
# for i in range(ngnr39):
#     domega = sol_domega[:,i]
#     plt.plot(dt, domega, "-", label="Omega Dot"+ str(i+1))
#
# plt.xlabel("Time")
# plt.ylabel("Omega Dot")
# plt.ylim(-5, 5)
# plt.legend();
# plt.title('Change Rates of Natural rotation frequencies')
# plt.savefig('domega_corcov.png')
#


check_times = 100
KK = 100  # repetition times
thres = np.array([0.2, 2])  # thres1 is for omega, thres2 is for omega_dot
# iid_disturbances = normaldisturbances(n=ngnr39, k=KK, sigma=sigma)
# # cor_disturbances = correlated_disturbances(delta=1, alpha=1.5, tau=1.5, n=ngnr39, L=redL39, k=KK)
# rates39 = model39.Simulation(KK, check_times, sigma, thres, t, nn, iid_disturbances)
# df39 = pd.DataFrame({'Node': node_list, 'RoCoF': rates39['vcheck_domega'],
#                      'AFV': rates39['vcheck_omega'], 'AV': rates39['vcheck_any']})
# # print(df39.to_latex(index=False))
# #
# df39.to_excel('tables\original_network_cor_100.xlsx')

# G3905 = multi_edge_removing(unG39, ([4, 5], [15, 16], [4, 14]))
# A3905, redL3905, redA3905 = kron_reduction(n39, ngnr39, G3905)
# model3905 = PowerNetworkSolver(theta0, omega0, redA3905, ngnr39, D, M, K, OMEGA)
# rates3905 = model3905.Simulation(KK, check_times, sigma, thres, t, nn, cor_disturbances)
# df3905 = pd.DataFrame({'Node': node_list, 'RoCoF': rates3905['vcheck_domega'],
#                      'AFV': rates3905['vcheck_omega'], 'AV': rates3905['vcheck_any']})

# df3905.to_excel('tuple_cor.xlsx')

# edge_list39 = list(unG39.edges)
# num_edges39 = getNumLines(getLines(subnetwork39))
# rate_list39 = np.zeros((num_edges39, 4))
# # print()
# j = 0
# for i in edge_list39:
#     temp_edge = np.array(i)
#     temp_G39 = edge_removing(unG39, temp_edge)
#     temp_A39, temp_redL39, temp_redA39 = kron_reduction(n39, ngnr39, temp_G39)
#     temp_model39 = PowerNetworkSolver(theta0, omega0, temp_redA39, ngnr39, D, M, K, OMEGA)
#     temp_rates39 = temp_model39.Simulation(KK, check_times, sigma, thres, t, nn, cor_disturbances)
#     temp_df39 = pd.DataFrame(
#         {'Node': node_list, 'RoCoF': temp_rates39['vcheck_domega'], 'AFV': temp_rates39['vcheck_omega'],
#          'AV': temp_rates39['vcheck_any']})
#     rate_list39[j, :] = temp_df39.mean(axis=0)
#     j += 1

KK = 1000
can_alpha = np.arange(0, 3.1, 0.2)
alpha_rates = np.zeros((len(can_alpha), 3))
for i in range(len(can_alpha)):
    temp_disturbances = correlated_disturbances(delta=1, alpha=can_alpha[i], tau=1.5, n=ngnr39, L=redL39, k=KK)
    # disturbances = normaldisturbances(n=ngnr39,k=1,sigma=sigma)
    # sol0 = np.pad(temp_disturbances[0], (ngnr39, 0), 'constant', constant_values=(0, 0))
    # single_sol = model39.solkuramoto(sol0, dt)
    # sol_domega = model39.getDotOmega(single_sol[:, :ngnr39], single_sol[:, ngnr39:], nn)
    temp_rates = model39.Simulation(KK, check_times, sigma, thres, t, nn, temp_disturbances)
    alpha_rates[i, 0] = np.mean(temp_rates['vcheck_omega'])
    alpha_rates[i, 1] = np.mean(temp_rates['vcheck_domega'])
    alpha_rates[i, 2] = np.mean(temp_rates['vcheck_any'])

can_delta = np.arange(0, 3.1, 0.2)
delta_rates = np.zeros((len(can_delta), 3))
for i in range(len(can_delta)):
    temp_disturbances = correlated_disturbances(delta=can_delta[i], alpha=1.5, tau=1.5, n=ngnr39, L=redL39, k=KK)
    # disturbances = normaldisturbances(n=ngnr39,k=1,sigma=sigma)
    # sol0 = np.pad(temp_disturbances[0], (ngnr39, 0), 'constant', constant_values=(0, 0))
    # single_sol = model39.solkuramoto(sol0, dt)
    # sol_domega = model39.getDotOmega(single_sol[:, :ngnr39], single_sol[:, ngnr39:], nn)
    temp_rates = model39.Simulation(KK, check_times, sigma, thres, t, nn, temp_disturbances)
    delta_rates[i, 0] = np.mean(temp_rates['vcheck_omega'])
    delta_rates[i, 1] = np.mean(temp_rates['vcheck_domega'])
    delta_rates[i, 2] = np.mean(temp_rates['vcheck_any'])
