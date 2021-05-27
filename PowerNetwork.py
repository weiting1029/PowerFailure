# !/usr/bin/env python
# coding: utf-8

# In[3]:
import copy
import matplotlib.pyplot as plt
import pypsa
from pypower.api import case39
import networkx as nx
import numpy as np
from numpy.linalg import inv
import sys


def networkTransform(case):
    ppc = case  # try case14() for a smaller network
    network = pypsa.Network()
    network.import_from_pypower_ppc(ppc)
    network.determine_network_topology()
    for sub_network in network.sub_networks.obj:
        pypsa.pf.calculate_B_H(sub_network)

    return network, sub_network


def getGenerators(network):
    generators = network.generators["bus"]
    return generators


def getBuses(sub_network):
    return sub_network.buses()[["Pd"]]


def getLines(sub_network):
    baseMVA = 100
    lines = sub_network.branches()
    lines.bus0 = lines.bus0.astype(int)  # transforms bus0 into an integer (it was a string)
    lines.bus1 = lines.bus1.astype(int)  # transforms bus1 into an integer (it was a string)
    lines.sort_values(by=['bus0', 'bus1'], inplace=True)  # sort the network lines based on buses they connect
    lines["susceptances"] = list(1 / (baseMVA * lines["x_pu_eff"]))  # susceptance calculated and added to the dataframe
    lines[["bus0", "bus1", "susceptances"]]  # printing only the relevant part of the dataframe

    return lines[["bus0", "bus1", "susceptances"]]


def getNumBuses(df_buses):
    return df_buses.shape[0]


def getNumLines(df_lines):
    return df_lines.shape[0]


def getUndGraph(df_buses, df_lines, network):
    generators = network.generators["bus"]
    n = df_buses.shape[0]
    node_list = np.arange(n) + 1
    weighted_edge_tuples = [tuple(x) for x in df_lines.to_numpy()]
    G = nx.DiGraph()
    G.add_nodes_from(node_list)
    # G.add_edges_from(edge_tuples)
    G.add_weighted_edges_from(weighted_edge_tuples)
    unG = nx.Graph(G)  # transform it into an undirected graph in network x
    ngnr = generators.shape[0]  # the number of generators

    return n, ngnr, unG


# %% define the kron reduction function
def kron_reduction(n, ngnr, unG):
    A = nx.to_numpy_matrix(unG)
    L = np.diag(np.squeeze(np.asarray(np.sum(A, axis=0)))) - A
    redL = L[(n - ngnr):n, (n - ngnr):n] - \
           np.transpose(L[:n - ngnr, n - ngnr:n]) @ inv(L[:n - ngnr, :n - ngnr]) @ L[:n - ngnr, n - ngnr:n]
    redA = -redL
    np.fill_diagonal(redA, 0)
    return A, redL, redA


# %% define the function of removing edges
def edge_removing(unG, re_edge):
    node1 = re_edge[0]
    node2 = re_edge[1]
    copyunG = copy.deepcopy(unG)
    copyunG.remove_edge(node1, node2)

    if len(list(nx.connected_components(copyunG))) == 1:
        return True, copyunG
    else:
        # print("The graph is not connected without edge ("+ str(node1) + ", " +
        #       str(node2) + ")")
        return False, unG


def multi_edge_removing(unG, args):
    new_unG = copy.deepcopy(unG)
    for edge in args:
        print(edge)
        connecting, new_unG = edge_removing(new_unG, edge)
    return new_unG

# define the function of moving more edges
