#!/usr/bin/env python
# coding: utf-8

# In[3]:
import copy
import matplotlib.pyplot as plt
import pypsa
#import data from pypower  
from pypower.api import case39
ppc = case39() #try case14() for a smaller network
network = pypsa.Network()
network.import_from_pypower_ppc(ppc)



network.determine_network_topology()
for sub_network in network.sub_networks.obj:
    pypsa.pf.calculate_B_H(sub_network)


# In[4]:


#printing only the relevant part of the buses dataframe
#namely the buses indices and the power they consume
sub_network.buses()[["Pd"]]



# In[5]:


#print the list of generators
generators = network.generators["bus"]
print(generators)


# In[7]:


baseMVA=100
lines = sub_network.branches()
lines.bus0 = lines.bus0.astype(int) #transforms bus0 into an integer (it was a string)
lines.bus1 = lines.bus1.astype(int) #transforms bus1 into an integer (it was a string)
lines.sort_values(by=['bus0','bus1'], inplace=True) #sort the network lines based on buses they connect
lines["susceptances"] = list(1/(baseMVA*lines["x_pu_eff"])) #susceptance calculated and added to the dataframe
lines[["bus0","bus1","susceptances"]] #printing only the relevant part of the dataframe


# # Transform the pypsa network object to the networkx network object

# In[11]:


import networkx as nx
import numpy as np
#to show the full matrix
import sys
np.set_printoptions(threshold=sys.maxsize)


# In[12]:


df_buses = sub_network.buses()[["Pd"]]
df_lines = lines[["bus0","bus1","susceptances"]]

n = df_buses.shape[0]#the number of nodes/buses
m = df_lines.shape[0]#the number of edges/lines

node_list =  np.arange(n)+1
#convert dataframe to tuples in order to be used to generate edges
weighted_edge_tuples = [tuple(x) for x in df_lines.to_numpy()]
# print(edge_tuples)


#generate the networkx object
G = nx.DiGraph()
G.add_nodes_from(node_list)
# G.add_edges_from(edge_tuples)
G.add_weighted_edges_from(weighted_edge_tuples)
unG = nx.Graph(G)#transform it into an undirected graph 
ngnr = generators.shape[0]#the number of generators 



#%% define the function of removing edges 
def edge_removing(unG,re_edge):
    """

    Parameters
    ----------
    unG : the undirected gragh object of Networkx 
        the original undirected graph 
    re_edge : an array with two elements
        the edge to be removed 

    Returns
    -------
    copyunG : undirected graph object of Networkx 
        without the edge to be removed 

    """
    node1 = re_edge[0]
    node2 = re_edge[1]
    copyunG = copy.deepcopy(unG) 
    copyunG.remove_edge(node1,node2)
    if(len(list(nx.connected_components(unG)))==1):
        return copyunG
    else:
        print("The graph is not connected without edge ("+ str(node1) + ", " + 
              str(node2) + ")") 
    
    


#%% define the kron reduction function 
from numpy.linalg import inv
def kron_reduction(unG, ngnr):
    A = nx.to_numpy_matrix(unG)
    L = np.diag(np.squeeze(np.asarray(np.sum(A,axis = 0))))-A
    redL = L[(n-ngnr):n,(n-ngnr):n] - np.transpose(L[:n-ngnr,n-ngnr:n])@inv(L[:n-ngnr,:n-ngnr])@L[:n-ngnr,n-ngnr:n]
    redA = -redL
    np.fill_diagonal(redA, 0)
    return A, redL, redA 



#%%
# # In[49]:

# #partition the Laplacian matrix
# # print(generators)

# A, redL,redA = kron_reduction(unG,ngnr)
# print(redL)
# print(redA)
# #visualize the adjacent matrix
# import scipy.sparse as sparse
# sparseA = sparse.csr_matrix(A)
# plt.spy(sparseA)














