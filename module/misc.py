def list_annotations():
    list_anno = ["sample","ZT", "Genotype", "run", "cell_type_final", "region_automap_name", "cell_class",]
    return list_anno

def cell_class():

    dict_temp = {
        'Tanycyte': 'Epithelial',
        'Ependymal': 'Epithelial',
        'Oligodendrocyte': 'Glial',
        'Endothelial': 'Epithelial',
        'Microglia': 'Glial',
        'Astro NT': 'Glial',
        'OPC': 'Glial',
        'Pericyte': 'Epithelial',
        'VLMC': 'Epithelial',
        'ABC' : 'Epithelial',
        'Astro TE': 'Glial',
        'Astro NT': 'Glial',
        'Oligodendrocyte': 'Glial',
        'Astro TE' : 'Glial',
        'Pericyte' :'Epithelial',
        'Ependymal' : 'Epithelial',           
        'Endothelial'   : 'Epithelial',
        'Tanycyte': 'Epithelial',      
        'OPC' : 'Glial',
        'VLMC': 'Epithelial',      
        'Choroid' :
        'Epithelial',
        }
    return dict_temp

def clock_genes_list():
    
    clock_genes = ['Arntl','Clock','Cry1','Cry2','Nr1d1',"Per1",'Per2','Per3','Rora','Rorb','Rorc',"Npas2"]
    
    return clock_genes