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

def genes_list(type_):
    
    dict_list = {
    'clock' : ['Arntl','Clock','Cry1','Cry2','Nr1d1',"Per1",'Per2','Per3','Rora','Rorb','Rorc',"Npas2"],

    'mitochondria' : ["Acacb","Acod1","Agxt2","Aldh1l2","Arg2","Asb9","Bbc3","Coq8a","Cps1","Cyp24a1","Ddit4","Decr1","Dhfr","Diablo",
                    'Elac2','Etnppl','Fam210b','Fmc1','Gcat','Gldc','Glud1','Gstk1','Hadhb','Hspa9','Idh2',"Isca1","Lrpprc","Mrm3","Mtif2","Mtus1",
                    "Nadk2","Nfu1",'Nubpl','Ogdh','Phykpl','Pmaip1','Pnkd','Ppp3r2','Ppp6c','Prorp','Rab32','Rsad1','Selenoo','Star','Tdh','Tufm',
                    "Vwa8","Acaa2","Glyat","Hibadh","Abcb8","Afg",'Bdh1','Chchd3','Cox7a1','Cox8b','Cpt2','Cyp11a1','Cyp11b2','Cyp27a1',
                    'Dhodh','Immt','Mcu','Mcub','Mcur1','Micu1','Mpc1','Mpc2',"Mtfp1","Ndufa12","Ndufaf6","Ndufs1","Ndufs4",'Ndufv2',"Opa1","Phb2",
                    "Pisd",'Sdha','Sdhb','Sdhc','Sfxn1','Slc25a20','Slc25a22','Slc25a28','Slc25a31','Slc25a51','Slc8b1',"Smdt1","Timm50","Ucp1",
                    "Ucp3","Aifm1","Coa7","Cpox","Cyct","Micu2","Htra2","Aldh2",'Sod2','Acss1','Ak4','Aldh4a1','Dglucy','Fdx1','Grpel1','Hadh',
                    'Hspd1','Maip1','Mmut','Ndufa9','Otc','Pde12','Pdk1','Pdk4','Pitrm1','Ppif','Prodh','Sardh','Sirt4','Twnk','Cs','Echs1','Ivd',
                    'Got2','Bcl2l1','Bcl2l2','Cyp27b1','Nrp1','Sfxn3','Timmdc1','Acsl1','Acsl4','Bak1','Bcl2','Cpt1c','Fundc1','Gdap1','Gk','Gpam',
                    'Gpat2','Hk2','Letmd1','Maoa','Maob','Mavs','Mfn1','Mfn2','Miga1','Miga2','Mtch2','Mtx1','Pgam5','Rhot2','Rmdn3','Synj2bp',
                    'Vdac2','Vps13a','Acsl6','Cpt1a','Hk1','Mtch1','Pink1','Bad','Sirt5','Polg','Bnip3'],



    }
    
    return dict_list[type_]