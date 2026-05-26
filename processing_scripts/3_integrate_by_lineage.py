# import rapids_singlecell as rsc
# import rmm
# import cupy as cp
# from rmm.allocators.cupy import rmm_cupy_allocator
import os
import scanpy as sc
import pandas as pd
import numpy as np
import scvi
import scanpy as sc
#from scvi.external import MRVI
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import gc
from matplotlib.colors import to_hex

# rmm.reinitialize(
#     managed_memory=True, # Allows oversubscription
#     pool_allocator=False, # default is False
#     devices=0, # GPU device IDs to register. By default registers only GPU 0.
# )
# cp.cuda.set_allocator(rmm_cupy_allocator)

import anndata as ad

import scipy.sparse as sp
import os


#from scvi.external import MRVI
import anndata as ad

import gc



BASE='/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/nobackup_output/subcluster2000_v7/'
adata2=sc.read_h5ad( '/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_new.h5ad.final.filtered.nohealthy')


hvg_number_target=1500
hvg_number=hvg_number_target
HVG_BATCH_KEY="sample_id"
MAX_EPOCHS=50
run_scanvi=True
run_mrvi=False

BATCH_SIZE=128
N_LATENTS=10
N_LAYERS=1
neighbor=20
mindist=0.1
SCANVI_LABELS_KEY="lvl5_annotation"  #"lvl5_annotation_new" #"lvl5_annotation2"
SCANVI_UNLABELLED= "KC"

import os



for x in adata2.obs["lvl1_new"].unique():
    if x in ["T"]:
        print(x)
        adata=adata2[adata2.obs["lvl1_new"]==x]
        x=x.replace("/","_")
        ADATA_PATH = BASE + f'adata_{x}.h5ad'
        # if os.path.exists(ADATA_PATH):
        #     print(f"Skipping {x}, already exists at {ADATA_PATH}")
        #     continue  # skip to next iteration
        adata.write(ADATA_PATH)

        # if adata.shape[0] < 10_000:
        #     neighbor=10
        # else:

        counts = adata.obs['sample_id'].value_counts()
        small = counts[counts < 50].index.tolist()
        mask = ~adata.obs['sample_id'].isin(small)
        adata = adata[mask].copy()
        adata.obs.sample_id.value_counts()
        adata.obs["sample_id"] = adata.obs["sample_id"].astype("category")
        adata.obs["sample_id"] = adata.obs["sample_id"].cat.remove_unused_categories()
        sc.pp.filter_genes(adata, min_cells=10)
        adata.shape



        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)



        def apply_qc_thresholds(adata, MIN_N_GENES, MAX_TOTAL_COUNT, MAX_PCT_MT, label, MIN_TOTAL_COUNT=0,):
            """
            Apply thresholds to generate QC column which says if passed all
            """
            ## Cell cycle gene list
            cc_genes_csv=pd.read_csv('/nfs/team298/ls34/csv_files/cc_genes.csv',  names=["gene_ids"], skiprows=1)
            cc_genes_csv = cc_genes_csv["gene_ids"]
            cc_genes_csv = list(cc_genes_csv)
            # Mark MT/ribo/Hb/cell cycle genes
            adata.var['mt'] = adata.var_names.str.startswith('MT-')  
            adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
            adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]")) 
            #adata.var["hb"] = adata.var_names.str.startswith(("HBA1", "HBA2", "HBB", "HBD","HBM", "HBZ", "HBG1", "HBG2", "HBQ1"))
            adata.var["cc_fetal"] = adata.var_names.isin(cc_genes_csv)

            # Calculate QC metrics
            sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False) #percent_top=[20],

            conditions = [
                (adata.obs['n_genes_by_counts'] < MIN_N_GENES),
                (adata.obs['pct_counts_mt'] > MAX_PCT_MT),
                (adata.obs['total_counts'] > MAX_TOTAL_COUNT),
                (adata.obs['total_counts'] < MIN_TOTAL_COUNT),
                (adata.obs['pct_counts_mt'] <= MAX_PCT_MT) & (adata.obs['n_genes_by_counts'] >= MIN_N_GENES) & 
                (adata.obs['total_counts'] <= MAX_TOTAL_COUNT)  & 
                (adata.obs['total_counts'] >= MIN_TOTAL_COUNT)
            ]
            label_suffix = label.split("_")[-1]
            print(label_suffix)
            pass_name = "Pass_" + label_suffix
            values = ['Low_nFeature', 'High_MT', 'High total count', 'Low total count', pass_name]
            assert len(conditions) == len(values), "Mismatch between conditions and values for np.select"
            assert all(isinstance(v, str) for v in values), "np.select values must all be strings"
            adata.obs[label] = np.select(conditions, values, default='Fail')
            adata.obs[label] = adata.obs[label].astype('category')

            print(adata.obs[label].value_counts())
        apply_qc_thresholds(adata, MIN_N_GENES=200, MAX_TOTAL_COUNT=300_000, MAX_PCT_MT=15,  MIN_TOTAL_COUNT=1000, label="QC_mid")


        hypoxia = ["VEGFA",
        "TF",
        "SLC2A1-AS1",
        "FOXN1",
        "VDAC1",
        "ASMT",
        "PLS3",
        "GPI",
        "DARS",
        "SNAPC1",
        "SEC61G",
        "GTF2IRD2B",
        "SAP30",
        "ZMYND8",
        "RSBN1",
        "BNIP3L",
        "GTF2IRD2",
        "STC2",
        "NARF",
        "HK2",
        "INHA",
        "PCF11",
        "CBWD3",
        "RAD51-AS1",
        "S100P",
        "HIF1A",
        ]

        additional_genes_to_exclude = [                             'JUND', 'HSPA1A', 'DNAJB1', 'EEF1A1', 'HSP90AA1', 'FTH1', 'FTL', 'HSPB1', 'XIST', 'VGLL3', "MEG3",
                                      "JUNB", "HSPA1B",  "FOSB", "HSP90AA1", "FOS", "DNAJB4", 'HSPA6', 'JUN', "NEAT1", "SOD2", "SOD3", "G0S2", "MYC"] 


        #original_hvg = str(hvg_number) 
        additional_genes_to_exclude = additional_genes_to_exclude + hypoxia

        mask_to_exclude = (
            adata.var.cc_fetal | 
            adata.var.hb | 
            adata.var.mt |
           # adata.var.mt2 |
            #adata.var.col |
            adata.var.ribo |
            adata.var.index.isin(additional_genes_to_exclude)
        )
        mask_to_include = ~mask_to_exclude
        adata  = adata[:, mask_to_include]
        adata.shape
        sc.pp.highly_variable_genes(adata,  
                                n_top_genes=hvg_number_target, 
                                subset=False,
                                batch_key=HVG_BATCH_KEY,
                                check_values=False,
                               )  
        var_genes_all = adata.var.highly_variable
        var_genes_batch = adata.var.highly_variable_nbatches > 2
        var_select = adata.var.highly_variable_nbatches >= 2
        var_genes = var_select.index[var_select]
        hvg_number = len(var_genes)
        print(f"selected {hvg_number} HVGs!")
        import gc
        gc.collect()
        label_dict1 = adata.var['highly_variable_nbatches'].to_dict()
        #adata2=sc.read_h5ad(adata_path)
        label_dict2 = adata.var['highly_variable'].to_dict()

        del(adata)
        import gc
        gc.collect()

        adata=sc.read_h5ad(ADATA_PATH)
        #adata = adata[adata.obs["combined_anno"]!="Nonspecific"]
        adata.var['highly_variable_nbatches'] = adata.var.index.map(label_dict1).fillna(np.nan)
        adata.var['highly_variable'] = adata.var.index.map(label_dict2).fillna(False)
        print(f"Added HVGs")

        best_HVG_BATCH_MINIMUM = None
        closest_hvg_number = None
        closest_difference = float('inf')
        for HVG_BATCH_MINIMUM in [1,2,3,4,5,6,7,10,20,30,40,50, 60,70, 90,100,110, 120,135, 150,160,180, 200,220,250,300,320,350,375,400,450]:
            var_genes_batch = adata.var.highly_variable_nbatches > HVG_BATCH_MINIMUM
            var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM
            var_genes = var_select.index[var_select]
            hvg_number = len(var_genes)

            difference = abs(hvg_number - hvg_number_target)

            if difference < closest_difference:
                closest_difference = difference
                closest_hvg_number = hvg_number
                best_HVG_BATCH_MINIMUM = HVG_BATCH_MINIMUM
        HVG_BATCH_MINIMUM=best_HVG_BATCH_MINIMUM
        hvg_number=closest_hvg_number
        CAT_COVS=[]
        CAT_COVS_TEMP = [x.replace("_", "").lower() for x in CAT_COVS] 
        collapsed_string = "_".join(CAT_COVS_TEMP)
        if len(CAT_COVS) == 0:
            model_details= "HVGNUMBER" + str(hvg_number) + "__MINBATCH" + str(HVG_BATCH_MINIMUM) + "__MAXEPOCHS" + str(MAX_EPOCHS) + "__BATCHKEY" + HVG_BATCH_KEY
        else:
            model_details= "HVGNUMBER" + str(hvg_number) + "__MINBATCH" + str(HVG_BATCH_MINIMUM) + "__MAXEPOCHS" + str(MAX_EPOCHS) + "__COVS" + collapsed_string
        print(f"selected {hvg_number} HVGs!")
        var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM
        adata = adata[:, var_select]
        print(f"{hvg_number} selected -> {adata.shape}, batch number {HVG_BATCH_MINIMUM}")



        adata=adata.copy()
        if run_scanvi==True:
            print("RUN SCANVI")
            def run_scvi(adata_hvg,  hvg_number , max_epochs, batch_size_vae, N_LATENT=10, N_LAYERS=1):
                DISPERSION = 'gene-batch'
                try:
                    details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
                except:
                    details="missing"
                adata_save_name = 'umap_' + details +"__1"
                print(adata_save_name)
                scvi.model.SCANVI.setup_anndata(adata_hvg, 
                                         batch_key=HVG_BATCH_KEY,
                                          labels_key=SCANVI_LABELS_KEY,
                                                unlabeled_category=SCANVI_UNLABELLED,
                                              categorical_covariate_keys=["tech"],

                                            #    layer="counts"
                                               )
                model = scvi.model.SCANVI(adata_hvg, 
                                dispersion=DISPERSION,
                                n_latent = N_LATENT, 
                                n_layers = N_LAYERS,
                               )
                model.train(accelerator ='gpu', 
                            max_epochs=max_epochs,             
                            early_stopping=True,
                           early_stopping_patience=5,
                           batch_size=batch_size_vae)
                print("model trained")
                latent = model.get_latent_representation() 

                try:
                    count=1
                    plt.subplots(figsize=(10, 10))
                    for key in model.history.keys():
                        plt.subplot(4,3,count)
                        plt.plot(model.history[key])
                        plt.title(key)
                        count+=1
                    plt.show()    
                except: 

                    print("Error with count")
                    try:
                        print(count)
                    except:
                        print("can't print count")
                return adata_hvg, model
        elif run_scanvi==False:
            if run_mrvi==False:
                print("RUN scvi")
                #sample_key = "sample_id"
                def run_scvi(adata_hvg,  hvg_number , max_epochs,  batch_size_vae, N_LATENT=10, N_LAYERS=1):
                    DISPERSION = 'gene-batch'
                    try:
                        details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
                    except:
                        details="missingdetails"
                    adata_save_name = 'umap_' + details +"__1"
                    print(adata_save_name)
                    scvi.model.SCVI.setup_anndata(adata_hvg,  
                                                    batch_key=HVG_BATCH_KEY, 
                                                  categorical_covariate_keys=["tech"],
                                             #     layer="counts"

                                                   )
                    model = scvi.model.SCVI(adata_hvg, 
                                   )
                    model.train(max_epochs=max_epochs,             
                                early_stopping=True,
                                accelerator='gpu',
                               early_stopping_patience=5,  
                               batch_size=batch_size_vae)
                    print("model trained")
                    return adata_hvg, model
            else:
                print("mrvi not set up")

        adata, model_test = run_scvi(adata, 
                             hvg_number=hvg_number, 
                                      max_epochs=MAX_EPOCHS, 
                                  batch_size_vae=BATCH_SIZE,
                                                   N_LATENT=N_LATENTS,
                                                  N_LAYERS=N_LAYERS )
        print("trained. now load adata")
        try:
            predictions = model_test.predict(adata)
        except:
            print("predictions didn't work")
        adata=sc.read_h5ad(ADATA_PATH)
        #¢adata = adata[adata.obs["combined_anno"]!="Nonspecific"]

        print(adata.shape)
        try:
            adata.obs["scanvi_predictions2"]=predictions
        except:
            print("fail")
        latent = model_test.get_latent_representation() 
        adata.obsm["X_scvi"] = latent

        neighbor_id = f"neighbor_{neighbor}"   
        print("start neighbours")
        sc.pp.neighbors(adata, use_rep = 'X_scvi', metric = "euclidean", n_neighbors=neighbor,key_added=neighbor_id)
        print("neighbours done")
        print("start umap")
        sc.tl.umap(adata, min_dist=mindist, neighbors_key = neighbor_id ) 



        adata.write(ADATA_PATH +".6000")
        print("saved", x)

        sc.tl.paga(adata, groups=SCANVI_LABELS_KEY,  neighbors_key=neighbor_id)
        sc.pl.paga(adata,
               color=SCANVI_LABELS_KEY,
               #neighbors_key=f"n_{N_NEIGHBORS}",
               add_pos=True,    # ← this computes adata_i.uns['paga']['pos']
               show=True)
        sc.tl.umap(adata, init_pos='paga',   min_dist=0.2,
                 neighbors_key = neighbor_id
                  )
        adata.write(ADATA_PATH + ".paga6000")
        del(adata)
        gc.collect()



