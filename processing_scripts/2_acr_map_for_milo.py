import scvi
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os

#adata_path='/nfs/team298/ls34/adult_skin/final_adatas/adata_scrna_freeze.h5ad'
#PATH = '/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_new.h5ad.final.filtered'
adata_path='/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_new.h5ad.final.filtered.scrna'
adata = sc.read_h5ad(adata_path)
#adata=adata[adata.obs["tech"]!="xenium"]


prep_hvgs=True
hvg_number=6000
HVG_NUMBER=6000

MAX_EPOCHS=50
HVG_BATCH_KEY = "sample_id"
HVG_BATCH_MINIMUM=60
run_scanvi=False

neighbor=20
mindist=0.1

N_LAYERS=2
N_LATENT=30

SCANVI_LABELS_KEY="lvl5_annotation"
SCANVI_UNLABELLED="New/unlabelled/excluded"

data_dir = '/nfs/team298/ls34/new_disease_atlas/milopy_final/'





# PART 1: PREP DATA FROM A GIVEN FILE PATH
adata_nonlesional=adata[adata.obs["atlas_status_simple2"]=="Nonlesional"]
adata_nonlesional.write("/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_A.h5ad")


# adata_lesionalall=adata[(adata.obs["atlas_status_simple"].isin(['Psoriasis_Lesional' ,
# 'Psoriasis_Nonlesional' ,
# 'Eczema_Lesional'     ,
# 'Eczema_Nonlesional' ,
#     'Eczema_Lesional_Fiskin', 'Eczema_Nonlesional_Fiskin'
#                                                               ]))|
#     (adata.obs["atlas_status"].isin(['Psoriasis_Nonlesional_Reynolds' ,
# 'Psoriasis_Nonlesional_Ma' ,
# 'Psoriasis_Nonlesional_Luc'     ,
# 'Eczema_Lesional_Other' ,
#     'Eczema_Lesional_Fiskin', 'Eczema_Nonlesional_Fiskin', 'Eczema_Lesional_Reynolds',
#                                      'Eczema_Nonlesional_Reynolds',
#                                                               ]))|         

 #                      ]
keep_status = [
    "Psoriasis_Nonlesional_Reynolds",
    "Psoriasis_Nonlesional_Ma",
    "Psoriasis_Nonlesional_Luc",
  #  "Eczema_Lesional_Other",
    "Eczema_Lesional_Fiskin",
    "Eczema_Nonlesional_Fiskin",
    "Eczema_Lesional_Reynolds",
    "Eczema_Nonlesional_Reynolds",
    'Psoriasis_Lesional_Ma','Psoriasis_Lesional_Luc',
'Psoriasis_Lesional_Reynolds',

]

adata_lesionalall = adata[adata.obs["atlas_status"].isin(keep_status)].copy()      

# adata_lesionalall.X=adata_lesionalall.layers["counts"].copy()
# del(adata_lesionalall.layers["counts"])
# print(adata_lesionalall.X[:5,:5].A)
adata_lesionalall.write("/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_PC.h5ad") 


# adata_lesionalall=adata[adata.obs["atlas_status_reynolds"].isin([#'Psoriasis_lesional' ,
# #'Psoriasis_nonlesional'  ,
# #'Eczema_lesional'     ,
# #'Eczema_nonlesional' ,
#                       "Psoriasis_lesional_other",                                          
#                                             "Psoriasis_nonlesional_other",                    
#                                                                 ])]
# adata_lesionalall=adata_lesionalall[adata_lesionalall.obs.dataset_id.str.startswith("Ma")]

# adata_lesionalall.write("/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_PC_Ma.h5ad") 


# Step 2: integrate atlas_A
if prep_hvgs==True:
    adata_path="/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_A.h5ad"
else:
    adata_path="/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_A_addedHVGs.h5ad"
#try:
adata=sc.read_h5ad(adata_path) 
#except:
#    print("need to prep HVGs")
#    prep_hvgs=True
#    adata_path="/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_A.h5ad"
#    adata=sc.read_h5ad(adata_path) 


if prep_hvgs:
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

        adata.obs[label] = np.select(conditions, values)
        adata.obs[label] = adata.obs[label].astype('category')

        print(adata.obs[label].value_counts())

    apply_qc_thresholds(adata, MIN_N_GENES=500, MAX_TOTAL_COUNT=300_000, MAX_PCT_MT=20,  MIN_TOTAL_COUNT=2000, label="QC_hi")



    hypoxia = ["VEGFA",
    #"#PTGS2",
    "TF",
    "SLC2A1-AS1",
    #"DK1",
    #"FOXN1",
    # "MMP9",
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
    #"FAM139A",
    "GTF2IRD2",
    #"C3ORF30",
    "STC2",
    "NARF",
    "HK2",
    "INHA",
    "PCF11",
    #"C9ORF30",
    "CBWD3",
    "RAD51-AS1",
    #"KIAA0195L",
    "#S100P",
    "HIF1A",
    ]

    additional_genes_to_exclude = [
        #"MMP14", 
        #"TNFAIP6", "ENO1",# "PDPN", "PTGES", "MMP2",
                                 'JUND', 'HSPA1A', 'DNAJB1', 'EEF1A1', 'HSP90AA1', 'FTH1', 'FTL', 'HSPB1', 'XIST', 'VGLL3', #"MEG3",
                                  "JUNB", "HSPA1B",  "FOSB", "HSP90AA1", "FOS", "DNAJB4", 'HSPA6', 'JUN', "NEAT1", "SOD2", "SOD3", "G0S2", "MYC"]  #HSPA1B FOSB 'DLK1', 'FABP5']


    original_hvg = str(HVG_NUMBER) + "select" + str(HVG_BATCH_MINIMUM)
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
    cells_per_sample = adata.obs['sample_id'].value_counts()    
    good_samples = cells_per_sample[cells_per_sample > 2_000].index
    adata = adata[adata.obs['sample_id'].isin(good_samples)].copy()
    print("Number of samples remaining: ")
    print(len(adata.obs['sample_id'].unique()))
    sc.pp.highly_variable_genes(adata,  
                                n_top_genes=HVG_NUMBER, 
                                subset=False,
                                batch_key=HVG_BATCH_KEY,
                                check_values=False,
                                #layer="normalized"
                               ) #1000 - 10_000    n_top_genes=3000,
    var_genes_all = adata.var.highly_variable
    var_genes_batch = adata.var.highly_variable_nbatches > HVG_BATCH_MINIMUM
    var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM
    var_genes = var_select.index[var_select]
    hvg_number = len(var_genes)
    print(f"selected {hvg_number} HVGs!")
    #adata.write(adata_path)

    ##

    

if prep_hvgs:
    label_dict = adata.var['highly_variable_nbatches'].to_dict()
    label_dict2 = adata.var['highly_variable'].to_dict()
    del(adata)
    import gc
    gc.collect()
    adata=sc.read_h5ad(adata_path)
    adata.var['highly_variable_nbatches'] = adata.var.index.map(label_dict).fillna(np.nan)
    adata.var['highly_variable'] = adata.var.index.map(label_dict2).fillna(False)
    #adata2.write("/nfs/team298/ls34/disease_atlas/milopy_final/adata_ACR_addedHVGs.h5ad")
    print(f"Added HVGs")
    counts = adata.obs[HVG_BATCH_KEY].value_counts()
    counts = counts[counts<50]
    adata=adata[~adata.obs[HVG_BATCH_KEY].isin(counts.index.to_list())]
    adata.shape
    adata.write("/nfs/team298/ls34/new_disease_atlas/milopy_final/adata_A_addedHVGs.h5ad")
    print("saved")








best_HVG_BATCH_MINIMUM = None
closest_hvg_number = None
closest_difference = float('inf')


for HVG_BATCH_MINIMUM in [1,2,3,4, 5, 10, 20, 25,30,35, 40,50, 55,60,65,70, 80,90,100, 120, 150, 200]:
    var_genes_batch = adata.var.highly_variable_nbatches > HVG_BATCH_MINIMUM
    var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM
    var_genes = var_select.index[var_select]
    hvg_number = len(var_genes)
    
    # Calculate the difference between the current hvg_number and 6000
    difference = abs(hvg_number - HVG_NUMBER)
    
    # Update the best HVG_BATCH_MINIMUM if this one is closer to 6000
    if difference < closest_difference:
        closest_difference = difference
        closest_hvg_number = hvg_number
        best_HVG_BATCH_MINIMUM = HVG_BATCH_MINIMUM
HVG_BATCH_MINIMUM=best_HVG_BATCH_MINIMUM
hvg_number=closest_hvg_number
var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM

CAT_COVS=[]
CAT_COVS_TEMP = [x.replace("_", "").lower() for x in CAT_COVS] 
collapsed_string = "_".join(CAT_COVS_TEMP)
try:
    if len(CAT_COVS) == 0:
        model_details= "HVGNUMBER" + str(hvg_number) + "__MINBATCH" + str(HVG_BATCH_MINIMUM) + "__MAXEPOCHS" + str(MAX_EPOCHS) + "__BATCHKEY" + HVG_BATCH_KEY
    else:
        model_details= "HVGNUMBER" + str(hvg_number) + "__MINBATCH" + str(HVG_BATCH_MINIMUM) + "__MAXEPOCHS" + str(MAX_EPOCHS) + "__COVS" + collapsed_string
except:
    model_details="model_unspecified"
print(f"{hvg_number} selected -> {adata.shape}. With {HVG_BATCH_MINIMUM} batches! ")
print(f"selected {hvg_number} HVGs (target {HVG_NUMBER})!")
adata = adata[:, var_select]
#hvg_number = adata.var.highly_variable.sum()
adata=adata.copy()


if run_scanvi==True:
    print("RUN SCANVI")
    def run_scvi(adata_hvg,  hvg_number , max_epochs, batch_size_vae,HVG_BATCH_KEY, N_LATENT=10, N_LAYERS=1):
        DISPERSION = 'gene'
        try:
            details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
        except:
            details="missing"
        adata_save_name = 'umap_' + details +"__1"
        print(adata_save_name)
        scvi.model.SCANVI.setup_anndata(adata_hvg, 
                                     #layer="counts",
                                   #categorical_covariate_keys=CATEGORICAL_COV,
                                  #continuous_covariate_keys=CONTINUOUS_COV,
                                 batch_key=HVG_BATCH_KEY,
                                  labels_key=SCANVI_LABELS_KEY,
                                        unlabeled_category=SCANVI_UNLABELLED
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
        # adata.obsm["X_scvi"] = latent
        # u_mde = scvi.model.utils.mde(latent)
        # adata.obsm["X_mde"] = u_mde
        # sc.pl.embedding(
        #     adata,
        #     basis="X_mde",
        #     color=["fine_annotation"],
        #     ncols=1,
        # )
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
    print("RUN scvi")
    sample_key = "sample_id"
    def run_scvi(adata_hvg,  hvg_number , max_epochs,  batch_size_vae, HVG_BATCH_KEY, N_LATENT=10, N_LAYERS=1):
        DISPERSION = 'gene-batch'
        try:
            details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
        except:
            details="missingdetails"
        adata_save_name = 'umap_' + details +"__1"
        print(adata_save_name)
        try:
            if adata_hvg.layers["counts"] is not None:
                scvi.model.SCVI.setup_anndata(adata_hvg, 
                                            layer="counts",
                                                batch_key=HVG_BATCH_KEY,
                                                #                                labels_key="broad_annotation",
                                                # unlabeled_category="New/unlabelled/excluded"
                                               )
            else:
                print(adata_hvg.X[:5,:5].A)
                scvi.model.SCVI.setup_anndata(adata_hvg, 
                                           # layer="counts",
                                                batch_key=HVG_BATCH_KEY,
                                                #                                labels_key="broad_annotation",
                                                # unlabeled_category="New/unlabelled/excluded"
                                               )
        except:
            scvi.model.SCVI.setup_anndata(adata_hvg, 
                                           # layer="counts",
                                                batch_key=HVG_BATCH_KEY,
                                                #                                labels_key="broad_annotation",
                                                # unlabeled_category="New/unlabelled/excluded"
                                               )
        print(adata_hvg.shape)
        print("start training")
        model = scvi.model.SCVI(adata_hvg, 
                        dispersion=DISPERSION,
                        n_latent = N_LATENT, 
                        n_layers = N_LAYERS,
                       )
        model.train(max_epochs=max_epochs,             
                    early_stopping=True,
                   # accelerator='gpu',
                   early_stopping_patience=3, #use_gpu =True, 
                  batch_size=batch_size_vae)
        print("model trained")
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
else:
    print(f"run scvi not defined. run_scanvi status is {run_scanvi}")


adata, model_test = run_scvi(adata, 
                 hvg_number=hvg_number, 
                          max_epochs=MAX_EPOCHS, 
                      #    HVG_BATCH_MINIMUM=HVG_BATCH_MINIMUM, 
                      batch_size_vae=512,
                          #CATEGORICAL_COV=CAT_COVS,
                         # CONTINUOUS_COV=[],
                          HVG_BATCH_KEY=HVG_BATCH_KEY,
                                       N_LATENT=N_LATENT,
                                      N_LAYERS=N_LAYERS )


#adata.write(f'/nfs/team298/ls34/disease_atlas/data/adata_inflammatlas_integrated{hvg_number}.h5ad', compression='gzip')
if run_scanvi==True:
    ref_model  = f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scanvi_{hvg_number}'
    model_test.save(ref_model,
                   save_anndata=True, 
                    overwrite=True) 


else:
    ref_model =f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}'
    model_test.save(ref_model,
                   save_anndata=True, 
                    overwrite=True) 
from datetime import datetime
now = datetime.now()
timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
print(f"Saved! Time: {timestamp} -> {ref_model}")




### PART 3: MAP ON REFERENCE
#'/nfs/team298/ls34/new_disease_atlas/model_scanvi5_MiloA_6518/'
data_dir = '/nfs/team298/ls34/new_disease_atlas/milopy_final/'

#ref_model= data_dir
#data_dir + f'model_A_scvi_{hvg_number}'
target_adata=sc.read_h5ad(data_dir+"adata_PC.h5ad")

try:
    if "counts" in target_adata.layers and target_adata.layers["counts"] is not None:
        target_adata.X=target_adata.layers["counts"].copy()
    else:
        print(target_adata.X[:5,:5].A)
        target_adata.layers["counts"] =   target_adata.X.copy()
except:
    print("ERROR WITH COUNTS CHECK")
    print(target_adata.X[:10,:10].A)


vae_q = scvi.model.SCVI.load_query_data(
             target_adata,
             ref_model,
             inplace_subset_query_vars=True)

vae_q.train(max_epochs=MAX_EPOCHS, plan_kwargs=dict(weight_decay=0.0))# use_gpu='cuda')
print("trained")
vae_q.save(data_dir + f'model_A_scvi_{hvg_number}_QueryREYNOLDS',
                   save_anndata=True, overwrite=True)
latent = vae_q.get_latent_representation() 

print(f'Final output: /nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS')

import scanpy as sc
#adata = sc.read_h5ad('/nfs/team298/ls34/disease_atlas/milo_adatas/model_A_scvi_2447_QueryREYNOLDS/adata.h5ad')
#adata=sc.read_h5ad("/nfs/team298/ls34/disease_atlas/milo_adatas/adata_PC_reynolds_send2cl.h5ad")
adata=sc.read_h5ad(data_dir + f'model_A_scvi_{hvg_number}_QueryREYNOLDS/adata.h5ad')
adata.obsm["X_scarches"] = latent
adata.layers["counts"]=adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


neighbor_id = "neighbor_"+str(neighbor)   
sc.pp.neighbors(adata, use_rep = 'X_scarches', metric = "euclidean", n_neighbors=neighbor)
print("neighbours done")

sc.tl.umap(adata, min_dist=mindist) 
print(f"UMAP done. Min dist {mindist}")
adata.write(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS/adata_clustered.h5ad')



# # part 4: ma version
# target_adata=sc.read_h5ad(data_dir+"adata_PC_Ma.h5ad")
# try:
#     if "counts" in target_adata.layers and target_adata.layers["counts"] is not None:
#         target_adata.X=target_adata.layers["counts"].copy()
#     else:
#         print(target_adata.X[:5,:5].A)
#         target_adata.layers["counts"] =   target_adata.X.copy()
# except:
#     print("ERROR WITH COUNTS CHECK")
#     print(target_adata.X[:10,:10].A)
# vae_q = scvi.model.SCVI.load_query_data(
#         target_adata,
#         ref_model,
#         inplace_subset_query_vars=True)

# vae_q.train(max_epochs=100, plan_kwargs=dict(weight_decay=0.0))# use_gpu='cuda')
# vae_q.save(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryMA',
#                    save_anndata=True, overwrite=True)
# latent = vae_q.get_latent_representation() 

# print("SUCCESS")
# import scanpy as sc
# #adata = sc.read_h5ad('/nfs/team298/ls34/disease_atlas/milo_adatas/model_A_scvi_2447_QueryREYNOLDS/adata.h5ad')
# #adata=sc.read_h5ad("/nfs/team298/ls34/disease_atlas/milo_adatas/adata_PC_reynolds_send2cl.h5ad")
# adata=sc.read_h5ad(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryMA/adata.h5ad')

# adata.obsm["X_scarches"] = latent
# adata.layers["counts"]=adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# neighbor=20
# neighbor_id = "neighbor_"+str(neighbor)   
# sc.pp.neighbors(adata, use_rep = 'X_scvi', metric = "euclidean", n_neighbors=neighbor,key_added=neighbor_id)
# print("neighbours done")
# mindist=0.3
# sc.tl.umap(adata, min_dist=mindist, neighbors_key =neighbor_id) 
# print(f"UMAP done. Min dist {mindist}")
# adata.write(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryMA/adata_clustered.h5ad')





BASE_MODEL=f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS/'


adata=sc.read_h5ad(BASE_MODEL + 'adata_clustered.h5ad')

#keep = ["Healthy", "Nonlesional atlas", "Eczema_lesional", "Eczema_nonlesional"]
adata=adata[(adata.obs["atlas_status"].str.startswith("Eczema_L"))|
            (adata.obs["atlas_status"].str.startswith("Eczema_Non"))
             ]
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata.write(BASE_MODEL + 'adata_annotated_eczemaonly.h5ad')


"""
To do: 
KEEP =  ['Lyme_disease/Erythema migrans', 'Vitiligo', 'irAE']
adata = adata[adata.obs["Patient_status"].isin(KEEP)]FDS


"""




adata=sc.read_h5ad(BASE_MODEL + 'adata_clustered.h5ad')

adata=adata[(adata.obs["atlas_status"].str.startswith("Psoriasis_L"))|
            (adata.obs["atlas_status"].str.startswith("Psoriasis_Non"))
             ]
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


adata.write(BASE_MODEL + 'adata_annotated_psoriasisonly.h5ad')

# adata=sc.read_h5ad(BASE_MODEL + 'adata_clustered.h5ad')
# keep = ["Nonlesional atlas", "Eczema_lesional", "Psoriasis_lesional"]
# adata=adata[adata.obs["atlas_status_reynolds"].isin(keep)]
# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)
# adata.write(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS/adata_annotated_both.h5ad')







adata=sc.read_h5ad(BASE_MODEL + 'adata_clustered.h5ad')
#keep = ["Healthy", "Nonlesional atlas", "Eczema_lesional", "Psoriasis_lesional"]
#adata=adata[adata.obs["atlas_status_reynolds"].isin(keep)]
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.write(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS/adata_annotated_both.h5ad')







# BASE_MODEL2=f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryMA/'
# os.listdir(BASE_MODEL2)

# adata=sc.read_h5ad(BASE_MODEL2 + 'adata_clustered.h5ad')
# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

# adata.write(f'/nfs/team298/ls34/new_disease_atlas/milopy_final/model_A_scvi_{hvg_number}_QueryREYNOLDS/adata_annotated_psoriasisonlyMA.h5ad')

# print(f"Base model: {BASE_MODEL}")
# print(f"Base model2: {BASE_MODEL2}")


