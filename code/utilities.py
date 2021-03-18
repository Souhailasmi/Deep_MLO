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
    
    
    
def plot_degree_distribution(G, savefig = False, figname = ''):
    degrees = {}
    for node in G.nodes():
        degree = G.degree(node)
        if degree not in degrees:
            degrees[degree] = 0
        degrees[degree] += 1
    sorted_degree = sorted(degrees.items())

    deg = [k for (k,v) in sorted_degree]
    cnt = [v for (k,v) in sorted_degree]

    fig, ax = plt.subplots(figsize = (20,10))
    plt.bar(deg, cnt, width=0.80, color='b')
    plt.title("Degree Distribution")
    plt.ylabel("Frequency")
    plt.xlabel("Degree")
    ax.set_xticks([d+0.05 for d in deg])
    ax.set_xticklabels(deg)  
    if savefig:
        fig.savefig('./plots/'+figname)
    
    
    
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

def get_quadratic_laplacian_forms(df,vertices,edges):
    n=len(df)
    df = min_max_scaler(df.iloc[:,1:])
    signals = df.columns
    component_graph = nx.Graph() # for a directed graph use nx.DiGraph()
    component_graph.add_nodes_from(range(len(vertices)))  # add multiple nodes at once
    component_graph.add_edges_from(edges)
    L = nx.laplacian_matrix(component_graph).A
    quad = {}
    for signal in signals:
        f = np.array(df[signal].fillna(0).values)
        inner = np.inner(f, L.dot(f))
        quad[signal] = inner
        
        
    quadtratics = pd.DataFrame(data = quad,index =[0])
        
    return  quadtratics 

###########################################################################################


def get_high_level_graph(high_level, thresh=None):
    """
    Get the high level delaunay triangulation graph of connected regions.
    Each region is representend by a unique node.
    ---
    Arguments :
        high_level : dataframe obtained with get_tils
        thresh : maximum distance threshold for edges
    Output :
        coord : array, node coordinates
        edges: set of edges below the threshold
    """
    coord = []
    for i, row in high_level.iterrows():
        coord.append([row.cell_x_position, row.cell_y_position])
    coord = np.array(coord)
    tri = Delaunay(coord)
    edges = set()
    dists = []
    for tr in tri.vertices:
        for i in range(3):
            edge_idx0 = tr[i]
            edge_idx1 = tr[(i+1)%3]
            if (edge_idx1, edge_idx0) in edges:
                continue  # already visited this edge from other side
            p0 = coord[edge_idx0]
            p1 = coord[edge_idx1]
            dists.append(np.linalg.norm(p1 - p0))
            if thresh is not None :
                if np.linalg.norm(p1 - p0) <  thresh:
                    edges.add((edge_idx0, edge_idx1))
            else :
                edges.add((edge_idx0, edge_idx1))
    return coord, edges


############################################################################################


def project_signal(graph,signal, plot_signal= True, color='r', savefig=False,figname=''):
    laplacian_matrix = nx.laplacian_matrix(component_graph).A
    eigen_vals, eigen_vects = np.linalg.eig(laplacian_matrix)
    idx = (-eigen_vals).argsort()[::-1]   
    eigenValues = eigen_vals[idx]
    eigenVectors = eigen_vects[:,idx]
    signal = signal
    P_transpose = eigenVectors.transpose()
    projection = np.dot(P_transpose,signal)
    if plot_signal :
        fig = plt.figure(figsize=(15,7))
        plt.plot(eigenValues,projection, color='r')
        plt.grid()
        plt.ylabel("Magnitude", size="large")
        plt.xlabel("Nodes", size="large")
        plt.title(figname, size="large")
        plt.show();
        if savefig:
            fig.savefig('./plots/signals/'+ figname+'.png')
    return projection  


##############################################################################################


def plot_relplot(data ,x,y,color='Reds',title='',savefig=False):
    fig = sns.relplot(data=data,kind="line", x=x, y=y, height=5, aspect=3,palette=color)  
    plt.xlabel(x)
    plt.ylabel(y)
    plt.grid()
    plt.title(title)
    plt.show()
    if savefig :
        fig.savefig('./plots/'+title+'.png')
        
        
###############################################################################################

def plot_density_per_component(data,metric,components=[1,2,3],phenotypes=['T','B'], savefig=False):
    data_copy = data.copy()  
    
    filter_phenotype = data_copy[data_copy.apply(lambda x: (x.phenotype in phenotypes) and (x.component in components), axis=1)]

    fig = sns.displot(filter_phenotype,x=metric, hue="phenotype",col='component', kind="kde")
    if savefig:
        fig.savefig('./plots/density_per_component'+METRIC+'.png')
        
        
##############################################################################################

def order_phenotypes(pheno_1, pheno_2):
    if pheno_1 > pheno_2:
        return str(pheno_1) + '-' + str(pheno_2)
    else:
        return str(pheno_2) + '-' + str(pheno_1)
    
    
################################################################################################  

def scatter_plot(data, x, y, title='',  savefig=False):
    data_copy = data.copy()
    
    plt.figure(figsize=(15,7))
    cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
    ax = sns.scatterplot(data= data_copy,x=x, y=y, palette=cmap )
    plt.grid()
    plt.title(title)    
    if savefig:
        fig.savefig('./plot/'+ title)
        
##################################################################################################

def get_subnodes(G , attribute, value : list):
    
    nodes_list = []
    
    for n in G.nodes():
        if G.nodes[n][attribute] in value :
            
            nodes_list.append(n)
            
    return nodes_list 

##################################################################################################

def visualize_graph(G, with_labels=False, k=None, alpha=0.1, color_attribute ='phenotype' ,node_shape='.',bipartite = False, savefig =False, figname =''):
    #nx.draw_spring(G, with_labels=with_labels, alpha = alpha)
    fig = plt.figure(figsize=(30,20))

    tumor_nodes = get_subnodes(G, 'phenotype',['tumor'])
    color_mapping = {'tumor':'darkred', 'stroma':'peru', 'T':'springgreen', 
                     'B':'turquoise','macrophages':'khaki', 'dendtritic':'orange', 'NK':'gray', 'MISSING' :'slategray'}
    if bipartite:
        pos = nx.bipartite_layout(G, tumor_nodes ,align ='vertical')
     
    else:
        pos = nx.spring_layout(G)
    
    ec = nx.draw_networkx_edges(G, pos, alpha=alpha)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes() ,node_color=[color_mapping[G.nodes[n][color_attribute]] for n in G.nodes()], 
                             alpha=0.8, node_shape = node_shape)
    #plt.colorbar(nc)
    plt.axis('off')
    plt.show()
    if savefig:
        fig.savefig('./plots/'+figname)
        
######################################################################################################


def plot_component(graph, MIN_ELEMENTS_PER_CLUSTER= 200,component_number = 0 ,remove_stroma =False, bipartite =False, savefig=False, figname=''):
    
    components = list(nx.connected_components(graph))
    components = filtering_components(components, MIN_ELEMENTS_PER_CLUSTER)
    print('The cell graph contains', len(components), 'connected components')
    components = sorted(components, key=len, reverse=False)
    subgraphs = [cell_graph.subgraph(c).copy() for c in components]
    main_component = subgraphs[component_number]
    if remove_stroma :
        nodes_list =[]
        for n in main_component.nodes():
            if main_component.nodes[n]["phenotype"] != 'stroma' :
                nodes_list.append(n) 
        main_component =main_component.subgraph(nodes_list)
        
   # print(nx.info(main_component))
    visualize_graph(main_component,bipartite=bipartite, savefig=savefig, figname=figname)
    
    
###########################################################################################################    

def min_max_scaler(data):
    df = data.copy()
    columns_name = df.columns
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled, columns = columns_name)
    
    return df


###########################################################################################################
        