import pandas as pd
import numpy as np
import pickle
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from google.cloud import bigquery
from scipy.spatial import Delaunay
from numpy import linalg
from sklearn import preprocessing
import plotly.graph_objects as go
import random
import os
import seaborn as sns




def graph_wrangling(graph, thrd = 50, only_border = True):
    
    graph_copy = graph.copy()
    graph_copy = graph_copy[graph_copy.distance < thrd]
    
    if only_border:
        
        boolean_condition = graph_copy.tissue_category_1 != graph_copy.tissue_category_2 
        graph_copy = graph_copy[boolean_condition]
        
    return graph_copy 

########################################################################################################

def filtering_components(list_components,MIN_NODES=1):
    filtered_components = []
    for component in list_components :
        if len(component)> MIN_NODES :
            filtered_components.append(component)
    return filtered_components 


########################################################################################################

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
    
    
    
    
########################################################################################################    
    
def plot_degree_distribution(G, PLOT_PATH,savefig = False, figname = '',no_show=False):
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
        fig.savefig(PLOT_PATH + figname)
    if no_show:
        plt.close()    
    
    
    
########################################################################################################    
    
    

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
    
    elif marker == 'SOX10p_CD56p':
        return 'DP1'
    
    elif marker == 'CD11cp_CD68p':
        return 'DP2'
    
    
    else :
        return 'MISSING'

########################################################################################################


def Union(lst1, lst2): 
    final_list = list(set(lst1) | set(lst2)) 
    return final_list 

########################################################################################################

def get_cells_from_edges(df):
    
    cell_list_1 = list(dict.fromkeys(df.cell_id_1))
    cell_list_2 = list(dict.fromkeys(df.cell_id_2))
    cells_union = Union(cell_list_1 , cell_list_2)
    
    return cells_union

########################################################################################################

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


########################################################################################################

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


########################################################################################################

def project_signal(graph,signal, PLOT_PATH,plot_signal= True, color='r', savefig=False,figname='',no_show=False):
    laplacian_matrix = nx.laplacian_matrix(graph).A
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
        plt.ylabel("Coefficient", size="large")
        plt.xlabel("eigenverctors", size="large")
        plt.title(figname, size="large")
        if savefig:
            fig.savefig(PLOT_PATH+ figname+'.png')
        if no_show:
            plt.close()    
            
    return projection  



########################################################################################################

def plot_density_per_component(data,metric,PLOT_PATH,components=[1,2,3],phenotypes=['T','B'], savefig=False,no_show=False):
    data_copy = data.copy()  
    
    filter_phenotype = data_copy[data_copy.apply(lambda x: (x.phenotype in phenotypes) and (x.component in components), axis=1)]

    fig = sns.displot(filter_phenotype,x=metric, hue="phenotype",col='component', kind="kde")
    if savefig:
        fig.savefig(PLOT_PATH+'density_per_component'+METRIC+'.png')
    if no_show:
        plt.close()    
        
########################################################################################################


def order_phenotypes(pheno_1, pheno_2):
    if pheno_1 > pheno_2:
        return str(pheno_1) + '-' + str(pheno_2)
    else:
        return str(pheno_2) + '-' + str(pheno_1)
    
########################################################################################################  

def scatter_plot(data, x, y,PLOT_PATH, title='',  savefig=False, no_show =False):
        
    data_copy = data.copy()
    
    plt.figure(figsize=(15,7))
    cmap = sns.cubehelix_palette(dark=.3, light=.8, as_cmap=True)
    ax = sns.scatterplot(data= data_copy,x=x, y=y, palette=cmap )
    plt.grid()
    plt.title(title)    
    if savefig:
        ax.savefig(PLOT_PATH+ title)
    if no_show:
        plt.close()
        
########################################################################################################

def visualize_graph(G,PLOT_PATH, with_labels=False, k=None, alpha=0.1, color_attribute ='phenotype' ,node_shape='.',bipartite = False, savefig   =False, figname ='',no_show=False):
    
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
        fig.savefig(PLOT_PATH + figname)
    if no_show:
        plt.close()

########################################################################################################

def plot_component(graph,PLOT_PATH, MIN_ELEMENTS_PER_CLUSTER= 200,component_number = 0 ,remove_stroma =False, bipartite =False, savefig=False, figname='',no_show=False):
    
    components = list(nx.connected_components(graph))
    components = filtering_components(components, MIN_ELEMENTS_PER_CLUSTER)
    print('The cell graph contains', len(components), 'connected components')
    components = sorted(components, key=len, reverse=False)
    subgraphs = [graph.subgraph(c).copy() for c in components]
    main_component = subgraphs[component_number]
    if remove_stroma :
        nodes_list =[]
        for n in main_component.nodes():
            if main_component.nodes[n]["phenotype"] != 'stroma' :
                nodes_list.append(n) 
        main_component =main_component.subgraph(nodes_list)
        
   # print(nx.info(main_component))
    visualize_graph(main_component,PLOT_PATH=PLOT_PATH,bipartite=bipartite, savefig=savefig, figname=figname,no_show = no_show)
    
    
########################################################################################################    

def min_max_scaler(data):
    
    df = data.copy()
    columns_name = df.columns
    x = df.values #returns a numpy array
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df = pd.DataFrame(x_scaled, columns = columns_name)
    
    return df

########################################################################################################

def G_cross_function(edges, min_radius=10, max_radius=50, radius_numbers=20, TILs=['T','B']):
    
    ## Defining radius quantiles
    ## The maximum radius here is 50 , since the maximum distance in all dataframes is 50
    ## The maximum radius should be changed once we change the theresholding
    
    vector_distances =list(np.linspace(min_radius, max_radius,radius_numbers))
    df = edges.copy()
    n = len(vector_distances)
    distribution = {}
    d_0 = vector_distances[0]
    
    ## This condition is based ont the fact that we are working with the border cell
    
    TILs_condition = (df.phenotype_1.apply(lambda x: x in TILs) & df.phenotype_2.apply(lambda x: x == 'tumor')) | (df.phenotype_1.apply(lambda x: x == 'tumor') & df.phenotype_2.apply(lambda x: x in TILs))
    df = df[TILs_condition]
    total_number = len(df)
    if total_number == 0:
        return pd.DataFrame({})
    else:
        
        df_0 = df[df.distance< d_0].copy()

        distribution[d_0] = len(df_0)/total_number
    
        for index in range(n-1):
            d_2 = vector_distances[index+1]
            condition = df.distance< d_2
            df_res = df[condition].copy()
            distribution[d_2] = len(df_res)/total_number
            
        distribution = pd.melt(pd.DataFrame(distribution, index = [0])) 
        
        return  distribution
    
########################################################################################################

def compute_AUC(df):
    if len(df) == 0:
        return 0
    else :
        n = len(df)
        x = df.variable.values
        f_x = df.value.values
        area = x[0]*f_x[0]
    
        for index in range(n-1):
        
            dx = x[index+1] - x[index]
            df = f_x[index] 
            area += df*dx
        
        return area   


#########################################################################################################

def plot_relplot(data ,x,y,PLOT_PATH,color='Reds',title='',savefig=False, no_show = False):
    fig = sns.relplot(data=data,kind="line", x=x, y=y, height=5, aspect=3,palette=color)  
    plt.xlabel(x)
    plt.ylabel(y)
    plt.grid()
    plt.title(title)
    if savefig :
        fig.savefig(PLOT_PATH+title+'.png')
    if no_show:
        plt.close()    
               
    
    
########################################################################################################


def visualize_component_graph(G,PLOT_PATH, with_labels=False, k=None, alpha=1.0, node_shape='.',savefig=False,figname='',no_show=False):
    #nx.draw_spring(G, with_labels=with_labels, alpha = alpha)
    fig = plt.figure(figsize=(30,20))
    pos = nx.spring_layout(G, k=k)
    if with_labels:
        lab = nx.draw_networkx_labels(G, pos, labels=dict([(n, n) for n in G.nodes()]))
    ec = nx.draw_networkx_edges(G, pos, alpha=alpha)
    nc = nx.draw_networkx_nodes(G, pos, nodelist=G.nodes(), node_color='g', node_shape=node_shape)
    plt.axis('off')
    if savefig:
        fig.savefig(PLOT_PATH+figname)
    if no_show:
        plt.close()    
        
        
########################################################################################################        
        
    
def get_subnodes(G, attribute, value):
    
    nodes_list = []
    
    for n in G.nodes():
        if G.nodes[n][attribute] in value :
            
            nodes_list.append(n)
            
    return nodes_list    

########################################################################################################

high_volume_patients = {'1C1':False,'1D1':True,'1E1':True,'1E2':True,'1E3':True,'1E4':True,'1E5':True,'1E6':True,'1E':True,'1F1':False,'1J1':True ,'1J2':True,'1K1':False}

########################################################################################################


def get_infiltrated_tils(nodes, TILs = ['T','B','NK'], on_border = True):
    
    Total_number_tils = nodes[nodes.phenotype.apply(lambda x: x in TILs)].phenotype.value_counts()

    total_T = Total_number_tils['T']
    total_B = Total_number_tils['B']
    total_NK = Total_number_tils['NK']
    
    condition_tissue = (nodes.tissue_category == 'tumor')
    condition_phenotype = (nodes.phenotype.apply(lambda x: x in TILs))
    border_condition = (nodes.on_border== on_border)
    infiltrated_tils_border = nodes[condition_tissue & condition_phenotype & border_condition]
    df_border = pd.DataFrame(infiltrated_tils_border.phenotype.value_counts()).reset_index().rename(columns={'index':'TIL','phenotype':'number_infiltration'})
    df_border_pivotted = pd.pivot_table(df_border,columns='TIL').reset_index()
    df_border_pivotted.columns.name =None 
    for pheno in TILs:
        if not(pheno in df_border_pivotted.columns):
            df_border_pivotted[pheno] = 0
            
    df_border_pivotted['T'] = df_border_pivotted['T']/ total_T
    df_border_pivotted['B'] = df_border_pivotted['B']/ total_B
    df_border_pivotted['NK'] =  df_border_pivotted['NK']/ total_NK
    
    return df_border_pivotted[TILs]


########################################################################################################


def group_and_pivot(data, columns_to_group:list, count: str,fill_na =False):
    
    data_grouped = data.groupby(columns_to_group).agg('count')[count].reset_index()
    data_pivot_table = pd.pivot_table(data_grouped, values=count, index=[columns_to_group[0]],
                    columns=[columns_to_group[1]], aggfunc=np.sum).reset_index()
    data_pivot_table.columns.name = None
    if fill_na :
        data_pivot_table = data_pivot_table.fillna(0).reset_index()
        
    return data_pivot_table 


########################################################################################################


def get_df_degrees(edges_data, cell_data, HIGH_VALUES=False):
    
    if HIGH_VALUES  :
    
        degree_table_1 = group_and_pivot(edges_data,['cell_id_1','tissue_category_2'],'cell_id_2',False).rename(columns={'cell_id_1' : 'cell_id'})
        degree_table_2 = group_and_pivot(edges_data,['cell_id_2','tissue_category_1'],'cell_id_1',False).rename(columns={'cell_id_2' : 'cell_id'})
        degree_table = degree_table_1.merge(degree_table_2, on='cell_id', how='outer')
        degree_table = degree_table.fillna(0)
        degree_table['missing_degree'] = 0
        degree_table = degree_table.rename( columns = {'tumor':'tumor_degree','stroma':'stroma_degree'})
        degree_table['total_degree'] = degree_table.apply(lambda x :   x.tumor_degree + x.stroma_degree , axis =1)
    
        return degree_table
    
    else :

        degree_table_1 = group_and_pivot(edges_data,['cell_id_1','tissue_category_2'],'cell_id_2',False).rename(columns={'cell_id_1' : 'cell_id'})
        degree_table_2 = group_and_pivot(edges_data,['cell_id_2','tissue_category_1'],'cell_id_1',False).rename(columns={'cell_id_2' : 'cell_id'})
        degree_table = degree_table_1.merge(degree_table_2, on='cell_id', how='outer')
        degree_table = degree_table.fillna(0)

        tissues = cell_data.tissue_category.unique()
        for tissue in tissues:
            condition = (tissue in degree_table_1.columns) and (tissue in degree_table_2.columns)
            if condition :
                print(tissue)
                column_name = tissue + '_degree'
                degree_table[column_name] = degree_table.apply(lambda x :   x[tissue+'_x'] + x[tissue+'_y'] , axis =1) 
                degree_table = degree_table.drop([tissue+'_x',tissue+'_y'],axis =1)
        
            else :
        
                degree_table = degree_table.rename(columns = {tissue : tissue+'_degree'})
            
            
        if 'missing_degree' in  degree_table.columns:   
            degree_table['total_degree'] = degree_table.apply(lambda x : x.tumor_degree + x.stroma_degree+ x.missing_degree , axis =1)
            
        else:
            degree_table['total_degree'] = degree_table.apply(lambda x : x.tumor_degree + x.stroma_degree , axis =1)
    
        return degree_table

    
########################################################################################################    


def cut_cercular_region(nodes, center=[0,0], diameter = 100):
    
    nodes_x_coordinates = nodes.cell_x_position.values
    nodes_y_coordinates = nodes.cell_y_position.values
    nodes_id = nodes.cell_id.values
    
    in_cercle_nodes = []
    points = []
    
    for cell_id,x,y in zip(nodes_id,nodes_x_coordinates,nodes_y_coordinates):
        
        if np.sqrt(np.power((x-center[0]),2) + np.power((y-center[1]),2) ) < diameter :
            
            in_cercle_nodes.append(cell_id)
            points.append([x,y])
            
    return in_cercle_nodes, points

########################################################################################################


def set_query_coordinates(patient_id, x_lower =0, x_upper = 1000, y_lower=0, y_upper=1000, thrd =0):

    query =  '''SELECT
    distance,
    cell_id_1,
    cell_id_2,
    tissue_category_1,
     tissue_category_2,
     phenotype_1,
  phenotype_2,
 FROM (
  SELECT
    table_1.cell_id_1,
    table_2.cell_id_2,
    table_1.x_1,
    table_1.y_1,
    table_2.x_2,
    table_2.y_2,
    tissue_category_1,
    tissue_category_2,
    phenotype_1,
    phenotype_2,


    SQRT(POW((x_1-x_2),2)+POWER((y_1-y_2),2)) AS distance
  FROM (
    SELECT
      cell_x_position AS x_1,
      cell_y_position AS y_1,
      cell_id AS cell_id_1,
      tissue_category as tissue_category_1,
      phenotype as phenotype_1
    FROM
      `advance-sonar-306410.deepmelo.DEEPMEL_{patient_id}_cell_seg_data`
      Where ((cell_x_position>{x_lower_1} and cell_x_position<{x_upper_1})and (cell_y_position>{y_lower_1} and cell_y_position<{y_upper_1})) ) AS table_1
  CROSS JOIN (
    SELECT
      cell_x_position AS x_2,
      cell_y_position AS y_2,
      cell_id AS cell_id_2,
      tissue_category as tissue_category_2,
      phenotype as phenotype_2,
 
    FROM
      `advance-sonar-306410.deepmelo.DEEPMEL_{patient_id}_cell_seg_data` 
      Where ((cell_x_position>{x_lower_2} and cell_x_position<{x_upper_2})and (cell_y_position>{y_lower_2} and cell_y_position<{y_upper_2}))) AS table_2 )

    where distance < {thrd} and (cell_id_1 < cell_id_2)'''.format(patient_id =patient_id,
                                                                  x_lower_1= x_lower,
                                                                  x_upper_1= x_upper,
                                                                  y_lower_1=y_lower,
                                                                  y_upper_1 = y_upper,
                                                                  x_lower_2= x_lower,
                                                                  x_upper_2=x_upper,
                                                                  y_lower_2=y_lower,
                                                                  y_upper_2 = y_upper,
                                                                  thrd = thrd)

    return query


########################################################################################################


def group_and_pivot(data, columns_to_group:list, count: str,fill_na =False):
    
    data_grouped = data.groupby(columns_to_group).agg('count')[count].reset_index()
    data_pivot_table = pd.pivot_table(data_grouped, values=count, index=[columns_to_group[0]],
                    columns=[columns_to_group[1]], aggfunc=np.sum).reset_index()
    data_pivot_table.columns.name = None
    if fill_na :
        data_pivot_table = data_pivot_table.fillna(0).reset_index()
        
    return data_pivot_table 


########################################################################################################


def get_average_axis_list(axis_data):
    
    phenotypes = list(axis_data.Phenotype.unique())
    pheno_map ={'macrophages': [4,3,4,6,5,5,7], 'NK':[3,3,4,6,6,5,7]} 
    
    for phenotype in phenotypes :
        pheno_map[phenotype] = axis_data[axis_data['Phenotype'] == phenotype].average_axis.values
        
    return pheno_map   
        
    
########################################################################################################    

def generate_axis(pheno_map, phenotype):
    
    return random.choice(pheno_map[phenotype])


########################################################################################################


def add_null_column(df, actual_columns):
    
    df_copy = df.copy()
    df_columns = df_copy.columns
    for actual in actual_columns:
        if not(actual in df_columns):
            df_copy[actual] = 0
    return df_copy 



########################################################################################################

def regions(tissue_category, border):
    
    if border==True :
        return 'border'
    else:
        return tissue_category
    
    
#####################################################################################################    