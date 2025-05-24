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
        'Choroid' : 'Epithelial',
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

    'astrocyte' : ["Acsbg1","Aqp4","Cdh20","Clmn","Gfap","Gli3","Id2","Mapk4","Ntsr2","Pde7b","Rfx4","Rorb","Slc39a12"], #
    'CA1-ProS' : ["Arhgap12","Fibcd1","Sipa1l3","Wfs1"],
# ["2010300C02Rik","Arhgef28","Bcl11b","Bhlhe22","Cabp7","Cpne4","Igfbp4","Necab2","Prdm8","Strip2","Syndig1"], #CA2
# ["Cpne6","Epha4","Hat1","Neurod6","Npy2r","Nrp2","Shisa6"], #CA3
# ["Cdh9","Orai2","Prox1","Rasl10a","Tanc1"], #DG
# ["Acvrl1","Adgrl4","Car4","Cd93","Cldn5","Cobll1","Emcn","Fgd5","Fn1","Kdr","Ly6a","Mecom","Nostrin","Paqr5","Pecam1","Pglyrp1","Slfn5","Sox17","Zfp366"], #Endothelial
# ["Arhgap25","Cd300c2","Cd53","Cd68","Ikzf1","Laptm5","Siglech","Sla","Spi1","Trem2"],#Microglia
# ['Gjc3','Gpr17','Opalin','Sema3d','Sema6a','Sox10','Zfp536'], #Oligodendrocytes
# ['Acta2','Ano1','Arhgap6','Carmn','Cspg4','Fos','Gucy1a1','Inpp4b','Nr2f2','Pip5k1b','Plekha2','Pln','Sncg','Sntb1'], # Pericytes
# ['Aldh1a2','Col1a1','Col6a1','Cyp1b1','Dcn','Fmod','Gjb2','Igf2','Pdgfra','Ror1',"Slc13a4","Spp1"], #VLMC
# ["Chat","Crh","Igf1",'Penk','Pthlh','Sorcs3','Thsd7a','Vip'], #Vip interneurons


    }
    
    return dict_list[type_]