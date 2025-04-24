### import necessary libraries
import anndata as ad
import scanpy as sc
from IPython.display import display
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
from module.misc import cell_class


### Automatic initial annotation
def automatic_initial_annotation(adata_spatial, cluster_col):
    cont_tab = pd.crosstab(adata_spatial.obs[cluster_col], adata_spatial.obs['mmc:subclass_name'], normalize="index")
    cont_tab_class = pd.crosstab(adata_spatial.obs[cluster_col], adata_spatial.obs['mmc:class_name'], normalize="index")
    max_col_dict = cont_tab.T.idxmax(axis=0).to_dict()
    max_col_dict_class = cont_tab_class.T.idxmax(axis=0).to_dict()
    adata_spatial.obs['cell_type_auto'] = adata_spatial.obs[cluster_col].map(max_col_dict)
    adata_spatial.obs['cell_class_auto'] = adata_spatial.obs[cluster_col].map(max_col_dict_class)

    all_cell_type = adata_spatial.obs['cell_type_auto'].unique()
    list_cell_nb = range(0, len(all_cell_type))
    mapping_dict = dict(zip(all_cell_type,list_cell_nb))
    adata_spatial.obs['cell_type_newnum_auto'] = adata_spatial.obs['cell_type_auto'].map(mapping_dict)

    all_class_type = adata_spatial.obs['cell_class_auto'].unique()
    list_cell_nb = range(0, len(all_class_type))
    mapping_dict = dict(zip(all_class_type,list_cell_nb))
    adata_spatial.obs['cell_class_newnum_auto'] = adata_spatial.obs['cell_class_auto'].map(mapping_dict)

    mapping_dict
    return adata_spatial


def auto_subclustering2(adata_to_sub, all_types = 'all', Clusters_to_use = 'cell_type_newnum_auto', resolution = 0.1):
    
    adata_filter = adata_to_sub
    
    bar = progressbar.ProgressBar(maxval=len(adata_filter.obs[Clusters_to_use].unique()), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    tracker_=0

    ### Automatic subclustering
    if Clusters_to_use == 'cell_type_newnum_auto':
        mapping_dict_all = dict(zip(adata_to_sub.obs['cell_id'], adata_to_sub.obs['cell_type_auto']))
    elif Clusters_to_use == 'cell_type_newnum_auto_sub':
        mapping_dict_all = dict(zip(adata_to_sub.obs['cell_id'], adata_to_sub.obs['cell_type_auto_sub']))
    
    clusters_numbers = Clusters_to_use
    
    if all_types == 'all':
        all_types = adata_to_sub.obs[clusters_numbers].unique() ### Subcluster all clusters
    current_ = 0
    for cluster_to_sub in (all_types):
        total_cluster = len(all_types)
        current_ = int(cluster_to_sub) + 1
        print(f'Subclustering of cluster {cluster_to_sub} ({current_} / {total_cluster})')
        adata_subcluster = adata_filter[adata_filter.obs[clusters_numbers] == cluster_to_sub]
        sc.pp.pca(adata_subcluster)
        sc.pp.neighbors(adata_subcluster)
        sc.tl.umap(adata_subcluster)
        sc.tl.leiden(adata_subcluster, resolution = 0.1)

        clustering_method = 'leiden'
        n = len(adata_subcluster.obs[clustering_method].unique())
        print(f'Results: {n} clusters')

        # adata_subcluster.obs['new_cluster'] = 0
        adata_subcluster.obs['new_cluster2'] = adata_subcluster.obs[clustering_method].astype("str")
        adata_subcluster.obs['new_cluster_class2'] = adata_subcluster.obs[clustering_method].astype("str")

        cont_tab = pd.crosstab(adata_subcluster.obs[clustering_method], adata_subcluster.obs['mmc:subclass_name'], normalize="index")
        max_col_dict = cont_tab.T.idxmax(axis=0).to_dict()

        adata_subcluster.obs['cell_type_auto'] = adata_subcluster.obs['new_cluster2'].map(max_col_dict)

        # Create a dictionary to map old values to new values
        mapping_dict = dict(zip(adata_subcluster.obs['cell_id'], adata_subcluster.obs['cell_type_auto']))

        mapping_dict_all.update(mapping_dict)
        
        tracker_ +=1
        bar.update(tracker_)

    # Use .map() function to rename cell contents in 'col1' based on mapping dictionary
    adata_to_sub.obs['cell_type_auto_sub'] = adata_to_sub.obs['cell_id'].map(mapping_dict_all)

    all_cell_type = adata_to_sub.obs['cell_type_auto_sub'].unique()
    list_cell_nb = range(0, len(all_cell_type))
    mapping_dict = dict(zip(all_cell_type,list_cell_nb))
    adata_to_sub.obs['cell_type_newnum_auto_sub'] = adata_to_sub.obs['cell_type_auto_sub'].map(mapping_dict)


########

def cluster_table(adata_to_use, Clusters_to_use = 'cell_type_newnum_auto_sub', sort_order = 'Cell Count', sort_ascend = False):

    adata_to_use.obs[Clusters_to_use] = adata_to_use.obs[Clusters_to_use].astype(int)

    bar = progressbar.ProgressBar(maxval=len(adata_to_use.obs[Clusters_to_use].unique()), \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    tracker_=0


    cont_tab_sub = pd.crosstab(adata_to_use.obs[Clusters_to_use],adata_to_use.obs['mmc:subclass_name'], normalize="index")
    cont_tab_sub = cont_tab_sub.loc[:, cont_tab_sub.sum(axis=0) > 0.05]

    cont_tab = pd.crosstab(adata_to_use.obs[Clusters_to_use], adata_to_use.obs['mmc:class_name'], normalize="index")
    cont_tab = cont_tab.loc[:, cont_tab.sum(axis=0) > 0.1] 

    all_cluster = adata_to_use.obs[Clusters_to_use].unique().astype(int)
    cluster_df = []
    max_corr = []
    celltype_ = []
    cellclass_ = []
    count_celltype = []
    per_celltype = []
    total_cells = len(adata_to_use)

    for all in all_cluster:
        # print(all)
        coor_temp = cont_tab_sub.T[all].sort_values(ascending = False)[0]
        max_corr.append(coor_temp)
        celltype_temp = cont_tab_sub.T[all].sort_values(ascending = False).index[0]
        celltype_.append(celltype_temp)
        cellclass_temp = cont_tab.T[all].sort_values(ascending = False).index[0]
        cellclass_.append(cellclass_temp)
        count_temp = adata_to_use.obs[adata_to_use.obs[Clusters_to_use] == all].groupby(Clusters_to_use)["cell_id"].count().values[0]
        count_celltype.append(count_temp)
        perc_temp = count_temp / total_cells * 100
        per_celltype.append(perc_temp)
        tracker_ +=1
        bar.update(tracker_)

    pd.set_option('display.max_rows', 250)

    d ={"Correlation":max_corr,"Celltype":celltype_,"Cell Class":cellclass_, "Cell Count":count_celltype, "Percentage": per_celltype}
    cluster_df = pd.DataFrame(data = d)

    if sort_order != 'None':
        cluster_df = cluster_df.sort_values(by=sort_order, ascending=sort_ascend)

    return cont_tab, cont_tab_sub, cluster_df


######
from module.misc import cell_class

def cell_class_annotation(adata):
    adata.obs['cell_class'] = 'Neuronal'

    dict_type = dict(zip(adata.obs['cell_type_final'],adata.obs['cell_class']))
    dict_temp = cell_class()
    dict_type.update(dict_temp)
    adata.obs['cell_class'] = adata.obs['cell_type_final'].map(dict_type)

    return adata

#####
import progressbar
from module.misc import clock_genes_list

def circascore_annot(adata, df):
    '''
    Create annotation based on the expression of clock genes in each cells
    '''
    bar = progressbar.ProgressBar(maxval=len(df), \
        widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    tracker_=0

    df['circascore']=0

    clock_genes = clock_genes_list()

    for n in range(0, len(df)):
        count_ = 0
        for g in clock_genes:
            count_ +=1
            if df[g][n] > 0:
                df['circascore'][n] += 1
        tracker_ +=1
        bar.update(tracker_)

    df['cell_id'] = df.index
    score_dict = dict(zip(df['cell_id'], df['circascore']))
    adata.obs['circascore'] = df['cell_id'].map(score_dict)

    return adata