from powerNetworkSolver import PowerNetworkSolver
from powerNetwork import networkTransform, getGenerators, getBuses, getLines, getNumBuses, getNumLines
from powerNetwork import getUndGraph, kron_reduction, edge_removing
from disturbancesGnr import normaldisturbances
from violationChecking import globalcheck
import seaborn as sns


sns.set_theme()
# importing the function for multi-edge removal
from powerNetwork import multi_edge_removing
from copy import deepcopy
import numpy as np
import pandas as pd


def greedy_algorithm(graph, n, ngnr, int_theta, int_omega, D, M, K, OMEGA, KK, check_times, sigma, thres, t, nn,
                     dtb_gnr, max_itr, type_rate):
    new_graph = deepcopy(graph)
    k = 0
    node_list = np.arange(ngnr) + 1
    while k <= max_itr:
        new_edge_list = list(new_graph.edges)
        new_num_edges = new_graph.number_of_edges()
        rate_list = np.zeros((new_num_edges,4))
        new_edge = np.array(new_edge_list[0])  # initialize the edge
        new_rate = 1  # initialize the rate
        j = 0
        for i in new_edge_list:
            temp_edge = np.array(i)
            temp_G = edge_removing(new_graph, temp_edge)
            temp_A, temp_redL, temp_redA = kron_reduction(n, ngnr, temp_G)
            temp_model39 = PowerNetworkSolver(int_theta, int_omega, temp_redA, ngnr, D, M, K, OMEGA)
            temp_rates39 = temp_model39.Simulation(KK, check_times, sigma, thres, t, nn, dtb_gnr)
            temp_df39 = pd.DataFrame(
                {'Node': node_list, 'RoCoF': temp_rates39['vcheck_omega'], 'AFV': temp_rates39['vcheck_theta'],
                 'AV': temp_rates39['vcheck_any']})
            rate_list[j, :] = temp_df39.mean(axis=0)
            if (rate_list[j, type_rate] < new_rate):
                new_edge = temp_edge
                prev_rate = new_rate
                new_rate = rate_list[j, type_rate]
            j += 1
        prev_graph = deepcopy(new_graph)
        print(new_edge) #print the edge to be removed
        new_graph = edge_removing(prev_graph, new_edge)
        k += 1
        print(str(new_graph.number_of_edges())+"  " + str(k))
    return new_graph

# show edge list
# # print
# edge_list39 = list(unG39.edges)
# num_edges39 = getNumLines(getLines(subnetwork39))
# rate_list39 = np.zeros((num_edges39,4))
# print()
# j=0
# for i in edge_list39:
#     temp_edge = np.array(i)
#     temp_G39 = edge_removing(unG39, temp_edge)
#     temp_A39, temp_redL39,temp_redA39 = kron_reduction(n39,ngnr39,temp_G39)
#     temp_model39 = PowerNetworkSolver(theta0,omega0,temp_redA39,ngnr39,D,M,K,OMEGA)
#     temp_rates39 = temp_model39.Simulation(KK, check_times, sigma, thres, t, nn, normaldisturbances)
#     temp_df39 = pd.DataFrame({'Node': node_list,'RoCoF': temp_rates39['vcheck_omega'] , 'AFV': temp_rates39['vcheck_theta'] ,'AV': temp_rates39['vcheck_any']})
#     rate_list39[j,:]=temp_df39.mean(axis = 0)
#     j+=1
