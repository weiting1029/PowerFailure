# importing the function for multi-edge removal
from copy import deepcopy

import numpy as np
import pandas as pd
import seaborn as sns
import time
from powerNetwork import kron_reduction, edge_removing
from powerNetworkSolver import PowerNetworkSolver
import multiprocessing as mp

sns.set_theme()


def greedy_algorithm(graph, n, ngnr, int_theta, int_omega, D, M, K, OMEGA, KK, check_times, thres, t, nn,
                     max_itr, type_rate, disturbances, pool):
    new_graph = deepcopy(graph)
    k = 0
    tol = 1e-5
    node_list = np.arange(ngnr) + 1
    new_rate = 1
    # prev_rate = 1
    continuation = True
    # start_time =0
    while k <= max_itr and continuation:
        new_edge_list = list(new_graph.edges)
        new_num_edges = new_graph.number_of_edges()
        rate_list = np.zeros((new_num_edges, 4))
        new_edge = np.array(new_edge_list[0])  # initialize the edge
        # new_rate = 1  # initialize the rate
        prev_rate = new_rate
        j = 0
        for i in new_edge_list:
            temp_edge = np.array(i)
            connecting, temp_G = edge_removing(new_graph, temp_edge)

            if connecting:
                temp_A, temp_redL, temp_redA = kron_reduction(n, ngnr, temp_G)
                temp_model39 = PowerNetworkSolver(int_theta, int_omega, temp_redA, temp_redL, ngnr, D, M, K, OMEGA)
                start_time = time.time()
                temp_list = temp_model39.parallelized_analytical_sml(check_times, thres, t, nn, disturbances, pool)
                # start_time = time.time()
                end_time = time.time()
                temp = end_time - start_time
                # print("total calculation time is time.time(): %f " % temp)
                df_total = pd.DataFrame(
                    {'Node': node_list, 'RoCoF': np.zeros(ngnr), 'AFV': np.zeros(ngnr), 'AV': np.zeros(ngnr)})
                for df in temp_list:
                    df_total['RoCoF'] = df_total['RoCoF'] + df['RoCoF']
                    df_total['AFV'] = df_total['AFV'] + df['AFV']
                    df_total['AV'] = df_total['AV'] + df['AV']

                df_total['RoCoF'] = df_total['RoCoF'] / KK
                df_total['AFV'] = df_total['AFV'] / KK
                df_total['AV'] = df_total['AV'] / KK

                # temp_rates39 = temp_model39.analytical_Simulation(check_times, thres, t, nn, disturbances)
                # temp_df39 = pd.DataFrame(
                #     {'Node': node_list, 'RoCoF': temp_rates39['vcheck_domega'], 'AFV': temp_rates39['vcheck_omega'],
                #      'AV': temp_rates39['vcheck_any']})
                #
                rate_list[j, :] = df_total.mean(axis=0)
                if rate_list[j, type_rate] < new_rate:
                    new_edge = temp_edge
                    # prev_rate = new_rate
                    new_rate = rate_list[j, type_rate]
                j += 1
            else:
                new_rate = prev_rate
                j += 1

        # prev_rate = new_rate
        if new_rate >= prev_rate:
            continuation = False  # end the outer loop

        prev_graph = deepcopy(new_graph)
        # print(new_edge)  # print the edge to be removed
        connecting2, new_graph = edge_removing(prev_graph, new_edge)
        k += 1
        print(str(new_graph.number_of_edges()) + " : " + str(k) + ", and the edge is " +
              str(new_edge))

    # pool.close()
    # pool.join()
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
