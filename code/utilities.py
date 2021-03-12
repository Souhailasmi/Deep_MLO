import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from google.cloud import bigquery




def graph_wrangling(graph, thrd = 50, only_border = True):
    
    graph_copy = graph.copy()
    graph_copy = graph_copy[graph_copy.distance < thrd]
    
    if only_border:
        
        boolean_condition = graph_copy.tissue_category_1 != graph_copy.tissue_category_2 
        graph_copy = graph_copy[boolean_condition]
        
    return graph_copy 

###################################################################################


def filtering_components(list_components,MIN_NODES=1):
    filtered_components = []
    for component in list_components :
        if len(component)> MIN_NODES :
            filtered_components.append(component)
    return filtered_components 


####################################################################################

def describe_graph(G):
    print(nx.info(G))
    if nx.is_connected(G):
        print("Avg. Shortest Path Length: %.4f" %nx.average_shortest_path_length(G))
        print("Diameter: %.4f" %nx.diameter(G)) # Longest shortest path
    else:
        print("Graph is not connected")
        print("Diameter and Avg shortest path length are not defined!")
    print("Sparsity: %.4f" %nx.density(G))  # #edges/#edges-complete-graph
    # #closed-triplets(3*#triangles)/#all-triplets
    print("Global clustering coefficient aka Transitivity: %.4f" %nx.transitivity(G))
    
    
    
######################################################################################    
    
    
    
def plot_degree_distribution(G):
    degrees = {}
    for node in G.nodes():
        degree = G.degree(node)
        if degree not in degrees:
            degrees[degree] = 0
        degrees[degree] += 1
    sorted_degree = sorted(degrees.items())

    deg = [k for (k,v) in sorted_degree]
    cnt = [v for (k,v) in sorted_degree]

    fig, ax = plt.subplots(figsize = (15,7))
    plt.bar(deg, cnt, width=0.80, color='b')
    plt.title("Degree Distribution")
    plt.ylabel("Frequency")
    plt.xlabel("Degree")
    ax.set_xticks([d+0.05 for d in deg])
    ax.set_xticklabels(deg)    
    
    
    
########################################################################################    
    
    

def get_phenotype (marker):
    
    if marker == 'SOX10p':
        return 'tumor'
    
    elif 'CD3p' in marker:
        return 'T'
    
    elif marker == 'DAPIp':
        return 'stroma'
    
    elif marker == 'CD20p':
        return 'B'
    
    elif marker == 'CD56p':
        return 'NK'
    
    elif marker == 'CD11Cp':
        return 'dendtritic'
    
    elif marker == 'CD68p':
        return 'macrophages'
    
    else :
        return 'MISSING'


###########################################################################################

def Union(lst1, lst2): 
    final_list = list(set(lst1) | set(lst2)) 
    return final_list 

#############################################################################################

def get_cells_from_edges(df):
    
    cell_list_1 = list(dict.fromkeys(df.cell_id_1))
    cell_list_2 = list(dict.fromkeys(df.cell_id_2))
    cells_union = Union(cell_list_1 , cell_list_2)
    
    return cells_union

############################################################################################