def get_most_probable_call_pciseq(cellData):
    # create dataframe with most probable cell types
    import pandas as pd

    append_data = pd.DataFrame(columns=['Cell_Num', 'X', 'Y', 'ClassName', 'Prob'])
    
    # create dataframe with most probable cell types
    for i, cell in enumerate(cellData.Cell_Num):
        cell_names = cellData.ClassName[i]
        cell_prob = cellData.Prob[i]
        try:
            max_prob = max(cell_prob)
        except TypeError: 
            max_prob = max([cell_prob])
        try: 
            index = [i for i, j in enumerate(cell_prob) if j == max_prob]
        except TypeError: 
            index = [i for i, j in enumerate([cell_prob]) if j == max_prob]
        cellname = cell_names[index[0]]
        X = cellData.X[i]
        Y = cellData.Y[i]
        data = [cell, X, Y, cellname, max_prob]
        append_data.loc[i] = data
    return append_data

def plot_most_probable_celltype(data, grouping = 'ClassName', color_column = 'color', 
                                invert_axis = False, legend = False, axis = True, 
                                title = 'Most probable celltype', size_of_dots = 20):
    import matplotlib.pyplot as plt 
    import matplotlib as mpl
    import pandas as pd
    
    plt.rcParams["figure.figsize"] = (20,20)
    data_group = data.groupby(grouping)
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling
    for name, group in data_group:
        ax.scatter(group.X, group.Y, s=size_of_dots,label=name, color = list(group[color_column]))
    #plt.xlim([-1000, 40000])
    #plt.ylim([-1000,45000])
    if invert_axis == True:
        plt.gca().invert_yaxis()
    if legend == True:
        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=10)
    if axis == False:
        plt.axis('off')
    plt.title(title + '\n' + 'Count: ' + str(data.shape[0]), size = 30)
    mpl.rcParams['figure.dpi']= 100
    plt.show()
    
def plot_most_probable_celltype_individually(data, grouping = 'ClassName',color_column = 'color', 
                                             invert_axis = False, legend = False, axis = True, 
                                             size_of_dots = 30):
    import matplotlib.pyplot as plt 
    import matplotlib as mpl
    import pandas as pd
    import numpy as np
    
    all_cells = np.unique(data[grouping])
    data_group = data.groupby(grouping)
    for cells in all_cells: 
        plt.rcParams["figure.figsize"] = (20,20)
        mpl.rcParams['figure.dpi']= 100
        count_graph = data_group.get_group(cells).shape[0]
        intermediate_df = data_group.get_group(cells)
        mean_probability = intermediate_df['Prob'].mean()
        plt.rcParams["figure.figsize"] = (15,12)
        plt.scatter(data.X, data.Y, s=size_of_dots, marker='o', alpha=.1, color = 'grey')
        plt.scatter((data_group.get_group(cells)).X, (data_group.get_group(cells)).Y, s=size_of_dots, color = list(intermediate_df[color_column])[0])
        plt.axis('off')
        plt.gca().invert_yaxis()
        plt.title(cells+ '\n'+ 'count: ' + str(count_graph)+'\n'+'percentage:' + str((round((count_graph/(data.shape[0]))*100, 2)))+'%' + '\n'+'mean_prob: '+str(round(mean_probability,3)),  fontsize=20)
        plt.show()
        
def from_polygon_to_labels(data,  polygon_file, output_prefix, x = 'x', y = 'y'):
    """
    The input files needed are the pciSeq output data and annotated polygons for the layers
    """
    from shapely.geometry import Polygon, Point,LineString,  mapping
    import json
    import geojson
    import numpy as np
    
    with open(polygon_file,'r') as r:
        shapejson = geojson.load(r)
    
    #last element in the geojson file is the column axis, so ignore it for now
    layer_annotations = shapejson
    cell_is_in_layer = {p["name"]:[Polygon(p["coordinates"][0]).intersects(Point(a))  for a  in data[[x,y]].values] for ii,p in enumerate(layer_annotations)}
    
    data["layer"] = "outside_VISp"
    for k in cell_is_in_layer.keys():
        data.loc[cell_is_in_layer[k],["layer"]] = k
    #data.to_csv(output_prefix + '_all_cells_with_layer_labels.csv', index = False )
    return data

def plot_cell_type_distribution_cortex(data,layer_color, dpi = 200, data_layer_column = 'layer', data_cluster_column = 'ClassName', data_color_column = 'color', title = 'Distribution of cell types in cortex layers'): 
    
    
    """The point of this function is for me to be able to plot layer statistics in an easy and convenient way. 
    
    INPUT: 
    - The input *DATA* should be pandas dataframes
    - The annotation table should include a color code and a cluster column with names that match the names in the input *DATA*
    """
    
    
    import pandas as pd
    import numpy as np
    import scipy.stats as scipy
    import matplotlib.pyplot as plt 
    
    all_layers = list(np.unique(data[data_layer_column]))
    #all_layers.remove('outside_VISp')
    cluster_to_color = dict(zip(data[data_cluster_column], data[data_color_column]))
    layer_color = layer_color.set_index('layer')
    list_of_all_clusters  = list(np.unique(data[data_cluster_column]))
    df1 = pd.DataFrame(index = list_of_all_clusters)

    fig, axs = plt.subplots(6,1, figsize=(10, 10), facecolor='w', edgecolor='k')
    plt.rcParams['figure.dpi'] = dpi

    plt.rc('figure', figsize=(15, 1.5))  
    group_layer = data.groupby(data_layer_column)
    group_cluster = data.groupby(data_cluster_column)
    
    for i, layers in enumerate(all_layers):
        intermediate_df = group_layer.get_group(layers)
        count_list = []
        cluster_list = []
        
        for j, cluster in enumerate(list(np.unique(intermediate_df[data_cluster_column]))):
            group_layer_cluster = intermediate_df.groupby(data_cluster_column)
            intermediate_intermediate_df = group_layer_cluster.get_group(cluster)
            intermediate_original = group_cluster.get_group(cluster)
            cluster_list.append(cluster)
            count_int = len(intermediate_intermediate_df)/len(intermediate_original)
            count_list.append(count_int)

        intermediate_dict = dict(zip(cluster_list, count_list))
        difference = (list(set(list_of_all_clusters) - set(cluster_list)))

        for k, not_there in enumerate(difference):
            intermediate_dict[not_there] = 0

        dataframe_counts = pd.DataFrame.from_dict(intermediate_dict, orient='index')
        dataframe_counts = dataframe_counts.sort_index()
        dataframe_counts['color'] = dataframe_counts.index.map(cluster_to_color)

        df1.insert(i,layers,dataframe_counts[0])

        axs[i].set_ylim([0, 1])
        axs[i].bar(dataframe_counts.index, dataframe_counts[0],color = dataframe_counts['color']) 
        axs[i].set_facecolor(list(layer_color.loc[layers])[0])
        fig.suptitle(title, fontsize=20)
        axs[i].set_ylabel(layers)

        if layers == 'L1_polygon' or layers == 'L2_polygon' or layers == 'L3_polygon' or layers == 'L4_polygon' or layers == 'L5_polygon':   
            axs[i].get_xaxis().set_ticks([])
            axs[i].get_yaxis().set_ticks([])
        else:
            axs[i].get_yaxis().set_ticks([])
            plt.xticks(dataframe_counts.index,list_of_all_clusters, rotation='vertical')
    return df1

def heatmap_layer_count(data, cmap_to_use = 'YlGnBu', title = 'Layer count'):
    import matplotlib.pyplot as plt
    import seaborn as sb
    data = data.dropna()
    x_axis = data.columns
    y_axis = data.index

    plt.rcParams["figure.figsize"] = (10,10)
    #sb.set(font_scale=2)
    ax = sb.heatmap(data, xticklabels=x_axis, yticklabels=y_axis,cmap=cmap_to_use)
    ax.set_title(title)
    plt.show()
    
def confusion_in_calling_normalized(cellData, most_probable):
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt 
    import seaborn as sb
    cellxgene_prob = pd.DataFrame(index=np.arange(len(cellData)), columns=np.arange(len(np.unique(list(most_probable.ClassName)))))
    cellxgene_prob.columns=np.unique(list(most_probable.ClassName))
    for i, cell in enumerate(cellData.Cell_Num):
        for j, gene in enumerate(np.unique(list(most_probable.ClassName))):
            gene_count_dict = dict(zip(cellData.ClassName[i], cellData.Prob[i]))
            try: 
                count = gene_count_dict[gene]
            except KeyError:
                count = 0
            cellxgene_prob.loc[cellxgene_prob.index[i], gene] = count
                                  
    cellxgene_prob = cellxgene_prob.fillna(0)
    cellxgene_prob.insert(0,'ClassName', most_probable.ClassName)
    
    cellxgene_prob_index = cellxgene_prob.set_index('ClassName')
    cellxgene_prob_index_mean = cellxgene_prob_index.groupby(by=cellxgene_prob_index.index, axis=0).mean()
    
    sum_probabilites = pd.DataFrame(index=np.arange(len(np.unique(list(most_probable.ClassName)))), columns=np.arange(len(np.unique(list(most_probable.ClassName)))))
    sum_probabilites.columns=np.unique(list(most_probable.ClassName))
    
    NEW_cellxgene_prob = pd.DataFrame(index=np.arange(len(np.unique(list(most_probable.ClassName)))), columns=np.arange(len(np.unique(list(most_probable.ClassName)))))
    NEW_cellxgene_prob.columns=np.unique(list(most_probable.ClassName))
    for i, celltype in enumerate(np.unique(cellxgene_prob_index.index)):
        values = list(cellxgene_prob_index.groupby(by=cellxgene_prob_index.index, axis=0).get_group(celltype).sum(axis = 0))
        values_normalized = np.array(values)/max(values)
        NEW_cellxgene_prob[celltype] = list(values_normalized)
    
    
    x_axis = NEW_cellxgene_prob.index
    y_axis = NEW_cellxgene_prob.columns
    #cellxgene_prob = np.log2(cellxgene_prob+1)
    plt.rcParams["figure.figsize"] = (50,25)
    sb.set(font_scale=1)
    plt.rcParams['figure.dpi'] = 100
    ax = sb.heatmap(NEW_cellxgene_prob, xticklabels=x_axis, yticklabels=y_axis,cmap="YlGnBu")
    return NEW_cellxgene_prob, cellxgene_prob_index

def genes_driving_classification(cellData, cellxgene, polygon):
    
    from shapely.geometry import Polygon, Point,LineString,  mapping
    import json
    import geojson
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sb
    
    all_genes = [] 
    for i, cell in enumerate(cellData.Genenames):
        all_genes.extend(cell)
    all_genes = np.unique(all_genes)
    all_genes = list(all_genes)

    all_cells = cellxgene.dropna()
    
    with open(polygon,'r') as r:
        shapejson = geojson.load(r)
    
    #last element in the geojson file is the column axis, so ignore it for now
    layer_annotations = shapejson
    cell_is_in_layer = {p["name"]:[Polygon(p["coordinates"][0]).intersects(Point(a))  for a  in cellxgene[["x","y"]].values] for ii,p in enumerate(layer_annotations)}
    
    cellxgene["layer"] = "outside_VISp"
    for k in cell_is_in_layer.keys():
        cellxgene.loc[cell_is_in_layer[k],["layer"]] = k
        
    #color_dataframe_2 = pd.DataFrame.from_dict(reference_cluster_color_dict,orient='index').sort_index()
    all_cells_no_nan = all_cells.dropna()
    all_cells_no_nan_index_only_genes = all_cells_no_nan.filter(all_genes).transpose().astype(float)
    all_cells_no_nan_index_only_genes_mean = all_cells_no_nan_index_only_genes.groupby(by=all_cells_no_nan_index_only_genes.columns, axis=1).mean()
    all_cells_no_nan_index_only_genes_tranposed_mean_log2_tranposed = np.log2(all_cells_no_nan_index_only_genes_mean+1).transpose()

    x_axis = all_cells_no_nan_index_only_genes_tranposed_mean_log2_tranposed.columns
    y_axis = all_cells_no_nan_index_only_genes_tranposed_mean_log2_tranposed.index
    plt.rcParams["figure.figsize"] = (50,50)
    sb.set(font_scale=0.5)
    ax = sb.clustermap(all_cells_no_nan_index_only_genes_tranposed_mean_log2_tranposed, xticklabels=x_axis, yticklabels=y_axis,cmap="YlGnBu") #metric="correlation"
    #ax.set_title('HCA09_1')
    plt.rcParams['figure.dpi'] = 500

    plt.show()
    
def add_color(data,annotation_file, column_to_map = 'ClassName',
                group_tag = 'ClassName', color_tag = 'colors'):
    to_map = dict(zip(annotation_file[group_tag], annotation_file[color_tag]))
    data['color'] = data[column_to_map].map(to_map)
    return data

def create_cellxgene_matrix(cellData, most_probable_call):
    
    import numpy as np
    import pandas as pd
    # create cellxgene matrix
    all_genes = [] 
    for i, cell in enumerate(cellData.Genenames):
        all_genes.extend(cell)
    all_genes = np.unique(all_genes)
    all_genes = list(all_genes)
    all_genes_column = list(all_genes)
    
    cellxgene = pd.DataFrame(index=np.arange(len(most_probable_call)), columns=np.arange(len(all_genes)))
    cellxgene.columns=all_genes
    coordinates_x = []
    coordinates_y = []
    highest_probability = []
    for i, cell in enumerate(cellData.Cell_Num):
        intermediate_x = cellData.X[i]
        intermediate_y = cellData.Y[i]
        coordinates_x.append(intermediate_x)
        coordinates_y.append(intermediate_y)

        cell_prob = cellData.Prob[i]
        max_prob = max(cell_prob)
        highest_probability.append(max_prob)

        for j, gene in enumerate(all_genes):
            gene_count_dict = dict(zip(cellData.Genenames[i], cellData.CellGeneCount[i]))
            try: 
                count = gene_count_dict[gene]
            except KeyError:
                count = 0
            cellxgene.loc[cellxgene.index[i], gene] = count
    cellxgene.insert(0, 'subclass', most_probable_call.ClassName)
    cellxgene['x'] = coordinates_x
    cellxgene['y'] = coordinates_y
    cellxgene['highest_probability'] = highest_probability
    cellxgene_no_zero = cellxgene[cellxgene['subclass'] != 'Zero']
    cellxgene_no_zero_withIndex = cellxgene_no_zero.set_index('subclass')
    return cellxgene_no_zero_withIndex

def prepare_single_cell_data_for_pciseq(single_cell_data, genes_in_columns, name_of_gene_column,annotation_file,
                                        column_name_of_cluster_to_map,column_name_of_desired_cluster,
                                        filter_genes, genes,log2_transform):
    
    import pandas as pd
    import numpy as np
    
    if filter_genes == True and genes_in_columns == True:
        single_cell_data_filtered = single_cell_data[single_cell_data[name_of_gene_column].isin(genes)]
    else: 
        print("Please transpose matrix")
    single_cell_data_filtered = single_cell_data_filtered.rename(columns={name_of_gene_column: 'GeneNames'})
    single_cell_data_filtered = single_cell_data_filtered.set_index('GeneNames')
    dictionary = dict(zip(annotation_file[column_name_of_cluster_to_map], annotation_file[column_name_of_desired_cluster]))
    single_cell_data_filtered.columns = single_cell_data_filtered.columns.map(dictionary)
   
    if log2_transform == True:
        single_cell_data_filtered = np.log2(single_cell_data_filtered+1)

    return single_cell_data_filtered

def plot_gene_expression(expression_data, gene_column = 'gene', x_column = 'x', y_column = 'y', 
                        size_of_dots = 1, title = 'Gene expression', axis = True, 
                        invert_x_axis = False, invert_y_axis = False, 
                        figure_dimesion = (20,20), figure_dpi = 100):

    import matplotlib as mpl 
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import random
    from matplotlib import colors as mcolors


    plt.rcParams["figure.figsize"] = (20,20)
    mpl.rcParams['figure.dpi']= figure_dpi

    expressiondata_groups = expression_data.groupby(gene_column)
    fig, ax = plt.subplots()
    ax.margins(0.05) # Optional, just adds 5% padding to the autoscaling

    markers = ['o', 'v', '^', '<', '>', '1', '2','3' ,'4'  ,'8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']

    colors = list((mcolors.CSS4_COLORS).values())
    for i, group in enumerate(expressiondata_groups.groups.keys()):
        number = (random.randint(0,len(markers)-1))
        ax.scatter(expressiondata_groups.get_group(group)[x_column], expressiondata_groups.get_group(group)[y_column], 
                    s=size_of_dots,c = colors[i], marker = markers[number],label=group) 
    if invert_y_axis == True:
        plt.gca().invert_yaxis()

    if invert_x_axis == True:
        plt.gca().invert_xaxis()

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=10)
    plt.axis('off')
    plt.title(title + '\n'+'(most probable celltype)' + '\n' + 'Count: ' + str(expression_data.shape[0]), size = 30)
    plt.show()

def plot_gene_expression_in_napari(expression_data, gene_column = 'gene', x_column = 'x', y_column = 'y', 
                        size_of_dots = 30, title = 'Gene expression', axis = True, 
                        invert_x_axis = False, invert_y_axis = False, 
                        figure_dimesion = (20,20), figure_dpi = 100):
    '''
    this funcion is a bit problematic to run but can be useful if one wants to interact with the data
    the problem is that napari needs some special commands
    gui qt5
    %gui qt
    '''

    import matplotlib as mpl 
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    import random
    from matplotlib import colors as mcolors
    import napari
    from IPython.core.display import display, HTML
    display(HTML("<style>.container { width:100% !important; }</style>"))
    from shapely.geometry import Polygon, Point,LineString,  mapping
    import matplotlib.pyplot as plt
    import seaborn as sb
    import geopandas
    napari.gui_qt()

    v = napari.Viewer()

    plt.rcParams["figure.figsize"] = (20,20)
    mpl.rcParams['figure.dpi']= figure_dpi

    expressiondata_groups = expression_data.groupby(gene_column)

    markers = ['o','+','-','x']
    colors = list((mcolors.CSS4_COLORS).values())

    for i, gene in enumerate(expressiondata_groups.groups.keys()):
        print(gene)
        number = (random.randint(0,len(markers)-1))
        spots_intermediate = expressiondata_groups.get_group(gene)
        spots_intermediate = spots_intermediate.reset_index()
        gene = v.add_points(spots_intermediate[[x_column,y_column]], symbol=markers[number],
                                     name = gene,edge_color= [0,0,0,0], size = size_of_dots,
                                     face_color= colors[i], blending = "translucent", opacity = 0.8)


def add_polygon_labels_to_genes(data, x, y, polygon_file):
    """
    The input files needed are the pciSeq output data and annotated polygons for the layers
    """
    from shapely.geometry import Polygon, Point,LineString,  mapping
    import json
    import geojson
    import numpy as np
    
    with open(polygon_file,'r') as r:
        shapejson = geojson.load(r)
    
    #last element in the geojson file is the column axis, so ignore it for now
    layer_annotations = shapejson
    cell_is_in_layer = {p["name"]:[Polygon(p["coordinates"][0]).intersects(Point(a))  for a  in data[[x,y]].values] for ii,p in enumerate(layer_annotations)}
    
    data["layer"] = "outside_VISp"
    for k in cell_is_in_layer.keys():
        data.loc[cell_is_in_layer[k],["layer"]] = k
    #data.to_csv(output_prefix + '_all_cells_with_layer_labels.csv', index = False )
    return data

def plot_gene_type_distribution_cortex(data, data_layer_column, gene_column, layer_color, title='Distribution of gene expression in cortex layers'): 
    
    
    """The point of this function is for me to be able to plot layer statistics in an easy and convenient way. 
    
    INPUT: 
    - The input *DATA* should be pandas dataframes
    - The annotation table should include a color code and a cluster column with names that match the names in the input *DATA*
    """
    
    
    import pandas as pd
    import numpy as np
    import scipy.stats as scipy
    import matplotlib.pyplot as plt 
    from matplotlib import colors as mcolors
    all_layers = list(np.unique(data[data_layer_column]))
    #all_layers.remove('outside_VISp')
    colors = list((mcolors.CSS4_COLORS).values())
    colors = colors[:len(np.unique(expressiondata[gene_column]))]
    
    cluster_to_color = dict(zip(np.unique(data['target']), colors))
    layer_color = layer_color.set_index('layer')
    list_of_all_clusters  = list(np.unique(data[gene_column]))
    df1 = pd.DataFrame(index = list_of_all_clusters)
    
    fig, axs = plt.subplots(6,1, figsize=(20, 10), facecolor='w', edgecolor='k')
    #plt.rcParams['figure.dpi'] = 200

    plt.rc('figure', figsize=(15, 1.5))  
    group_layer = data.groupby(data_layer_column)
    group_cluster = data.groupby(gene_column)
    
    for i, layers in enumerate(all_layers):
        intermediate_df = group_layer.get_group(layers)
        count_list = []
        cluster_list = []
        
        for j, cluster in enumerate(list(np.unique(intermediate_df[gene_column]))):
            group_layer_cluster = intermediate_df.groupby(gene_column)
            intermediate_intermediate_df = group_layer_cluster.get_group(cluster)
            intermediate_original = group_cluster.get_group(cluster)
            cluster_list.append(cluster)
            count_int = len(intermediate_intermediate_df)/len(intermediate_original)
            count_list.append(count_int)

        intermediate_dict = dict(zip(cluster_list, count_list))
        difference = (list(set(list_of_all_clusters) - set(cluster_list)))

        for k, not_there in enumerate(difference):
            intermediate_dict[not_there] = 0

        dataframe_counts = pd.DataFrame.from_dict(intermediate_dict, orient='index')
        dataframe_counts = dataframe_counts.sort_index()
        dataframe_counts['color'] = dataframe_counts.index.map(cluster_to_color)
        
        df1.insert(i,layers,dataframe_counts[0])

        axs[i].set_ylim([0, 1])
        axs[i].bar(dataframe_counts.index, dataframe_counts[0],color = dataframe_counts['color']) 
        axs[i].set_facecolor(list(layer_color.loc[layers])[0])
        fig.suptitle(title, fontsize=20)
        axs[i].set_ylabel(layers)

        if layers == 'L1_polygon' or layers == 'L2_polygon' or layers == 'L3_polygon' or layers == 'L4_polygon' or layers == 'L5_polygon':   
            axs[i].get_xaxis().set_ticks([])
            axs[i].get_yaxis().set_ticks([])
        else:
            axs[i].get_yaxis().set_ticks([])
            plt.xticks(dataframe_counts.index,list_of_all_clusters, rotation='vertical')
    return df1

def get_cellxgene_matrix(geneData):
    
    import numpy as np
    import pandas as pd

    all_genes = np.unique(geneData['Gene'])
    all_cells = np.unique(geneData['neighbour'])

    cellxgene = pd.DataFrame(index=np.arange(len(all_cells)), columns=np.arange(len(all_genes)))
    cellxgene.columns=all_genes
    cellxgene.index=all_cells
    geneData_group = geneData.groupby('neighbour')

    for i, cell in enumerate(all_cells):

        dataframe_counts = pd.DataFrame(geneData_group.get_group(cell)['Gene'].value_counts())
        dict_counts = dict(geneData_group.get_group(cell)['Gene'].value_counts())
        list_of_genes = list(dataframe_counts.index)

        for j, gene in enumerate(list_of_genes):
            count = dict_counts[gene]
            cellxgene.loc[cellxgene.index[i], gene] = count
    cellxgene = cellxgene.fillna(0)
    return cellxgene

def get_probability_matrix(cellData):
    import numpy as np 
    import pandas as pd 
    flat_list = [item for sublist in list(cellData['ClassName']) for item in sublist]
    celltypes = list(np.unique(flat_list))
    cell_index = np.unique(cellData['Cell_Num'])
    probability_matrix = pd.DataFrame(columns = np.arange(len(celltypes)), index = np.arange(len(cell_index)))
    probability_matrix.columns = celltypes
    probability_matrix.index = cell_index
    for i, cell in enumerate(np.unique(cellData['Cell_Num'])): 
        current_dictionary = dict(zip(cellData.loc[i]['ClassName'], cellData.loc[i]['Prob']))
        for j, key in enumerate(list(current_dictionary.keys())):
            probability_matrix.loc[i, key] = current_dictionary[key]
    probability_matrix = probability_matrix.fillna(0)
    probability_matrix_tranposed = probability_matrix.transpose()
    max_celltype = probability_matrix_tranposed.idxmax(axis = 0)
    return probability_matrix, max_celltype