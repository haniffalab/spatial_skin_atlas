#from scvi.external import MRVI
import scvi
import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
from matplotlib.colors import to_hex
from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering
from scipy.spatial.distance import squareform


N_LATENT=10
N_LAYERS=1
NEIGHBOR=20

MIN_DIST=0.1
HVG_BATCH_KEY = "info_id6"
HVG_NUMBER = 1500
HVG_BATCH_MINIMUM=80
MAX_EPOCHS=10

DATASET= "final_xenium_paper3"


run_scanvi=False
SCANVI_LABELS_KEY="lvl5_annotation"
SCANVI_UNLABELLED="New/unlabelled/excluded"
run_mrvi=False
prep_HVGS_done = False
MERGE_DATASETS=False
INT_PREPARED=False

# adata_path='/nfs/team298/ls34/adult_skin/final_adatas/adata_xenium_freeze_plus3d.h5ad.september'

# adata_i=sc.read_h5ad(adata_path)
# adata_i.obs["lvl4_annotation"]=adata_i.obs["lvl4_annotation"].fillna(adata_i.obs["scanvi_predictions"])
# adata_new=sc.read_h5ad('/nfs/team298/beacon_ls34/loaded_adatas/adata_SEPT25.h5ad')
# adata=ad.concat([adata_i, adata_new], label="fig8_mapping", keys=["original", "sept25"], join='outer')
#adata=sc.read_h5ad('/nfs/team298/ls34/adult_skin/final_adatas/adata_xeniumonly.h5ad.final.filtered.withvalidation.clean')
adata_path = '/nfs/team298/ls34/adult_skin/final_adatas/adata_xeniumonly.h5ad.final.filtered.withvalidation.clean'




if INT_PREPARED:

    adata=sc.read_h5ad(adata_path+".int")
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
    print(f"selected {hvg_number} HVGs! Target {HVG_NUMBER}")
    ###

    adata2=sc.read_h5ad(adata_path)
    label_dict = adata.var['highly_variable_nbatches'].to_dict()
    adata2.var['highly_variable_nbatches'] = adata2.var.index.map(label_dict).fillna(np.nan)
    label_dict = adata.var['highly_variable'].to_dict()
    adata2.var['highly_variable'] = adata2.var.index.map(label_dict).fillna(False)
    adata2.write(adata_path + ".HVGs")
    print(f"Added HVGs")

    adata=adata2.copy()
    #adata2=0
    import gc
    del adata2  # Delete the variable reference
    gc.collect()  # Force garbage collection

else:
    if prep_HVGS_done == False:
        adata=sc.read_h5ad(adata_path)
        print("adata loaded")
    # #             print(adata.shape)
    # #             adata = adata[adata.X.sum(axis=1) > 0]  # Remove empty cells
    # #             adata = adata[:, adata.X.sum(axis=0) > 0]  # Remove empty genes
    # #             print(adata.shape)
    # #             # Compute total counts per sample
    # #             sample_counts = adata.obs["sample_id"].value_counts()

    # #             # Keep only samples with at least 50 counts
    # #             valid_samples = sample_counts[sample_counts >= 50].index

    # #             # Filter adata to keep only valid samples
    # #             adata = adata[adata.obs["sample_id"].isin(valid_samples), :].copy()
                #adata.write(f'/nfs/team298/ls34/new_disease_atlas/adata_{DATASET}_and_Atlas.h5ad')
            #adata=sc.read_h5ad(f'/nfs/team298/ls34/new_disease_atlas/adata_{DATASET}_and_Atlas.h5ad.int')






    #     #adata_path = '/nfs/team298/ls34/new_disease_atlas/adata_all_scvi5.countsonly_notstar.h5ad'
    #     #'/nfs/team298/ls34/new_disease_atlas/adata_all_scvi5.h5ad.countsonly' - with Tstar

    #     # ADATA_PATH_FOR_LABELS = "/nfs/team298/ls34/disease_atlas/mrvi/adata_inflamm_scanvi6_HVG8024_EPOCHS10.h5ad.v2.fordeconv.noHS.reintegratedSCANVI.updated10"
    #     # adata_labelled = sc.read_h5ad(ADATA_PATH_FOR_LABELS)
    #     # label_dict = adata_labelled.obs["lvl5_annotation_new3"].to_dict()
    #     # import pickle

    #     # Now `labels_dict` contains the loaded dictionary
    #     # save_path = "/nfs/team298/ls34/labels_lvl5_annotation.pkl"
    #     # with open(save_path, 'wb') as file:
    #     #    pickle.dump(label_dict, file)
    #     # adata.obs['lvl5_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")

    #     # adata.obs['lvl3_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")

    #     # adata.obs['lvl5_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")
    #     # adata.layers["counts"]=adata.X.copy()
    #     # adata.write(adata_path)
    #     # print(f"updated labels and saved to {adata_path}")



    #     # print("save adata")

    #     # adata_i = adata[(adata.obs["leiden_res0.2"]=="22")|
    #     #                (adata.obs["leiden_res0.2"]=="23")|
    #     #                 (adata.obs["leiden_res0.2"]=="2")
    #     #                ]
    #     # adata_i.write(base_dir + 'adata_scvi4_removedjunk_LC_MigDC_only.h5ad')
    #     # adata_i = adata[adata.obs["leiden_res0.2"]=="7"]
    #     # adata_i.write(base_dir + 'adata_scvi4_Schwannrelated.h5ad')

    #     try:
    #         adata.X = adata.layers["counts"].copy()
    #     except:
    #         adata.layers["counts"] = adata.X.copy()
    #         print(adata.X[:6,:6].A)
    #         print("set adata.X to be above as no counts layer")
    #         # del(adata.layers["counts"])

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)


        hypoxia = ["VEGFA",
        #"#PTGS2",
        "TF",
        "SLC2A1-AS1",
        #"DK1",
        "FOXN1",
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
        "S100P",
        "HIF1A",
        ]
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

        additional_genes_to_exclude = [
            #"MMP14", 
            #"TNFAIP6", "ENO1",# "PDPN", "PTGES", "MMP2",
                                     'JUND', 'HSPA1A', 'DNAJB1', 'EEF1A1', 'HSP90AA1', 'FTH1', 'FTL', 'HSPB1', 'XIST', 'VGLL3', "MEG3",
                                      "JUNB", "HSPA1B",  "FOSB", "HSP90AA1", "FOS", "DNAJB4", 'HSPA6', 'JUN', "NEAT1", "SOD2", "SOD3", "G0S2", "MYC",
                                "MKIF7", "TOP2A",
        ]  #HSPA1B FOSB 'DLK1', 'FABP5']


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

        adata.write(adata_path + ".int")
        import gc

        del adata  # Delete the variable reference
        gc.collect()  # Force garbage collection
        adata=sc.read_h5ad(adata_path + ".int")
        cells_per_sample = adata.obs[HVG_BATCH_KEY].value_counts()       
        good_samples = cells_per_sample[cells_per_sample > 2_000].index

        adata = adata[adata.obs[HVG_BATCH_KEY].isin(good_samples)].copy()
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
        print(f"selected {hvg_number} HVGs! Target {HVG_NUMBER}")
        ###

        adata2=sc.read_h5ad(adata_path)
        label_dict = adata.var['highly_variable_nbatches'].to_dict()
        adata2.var['highly_variable_nbatches'] = adata2.var.index.map(label_dict).fillna(np.nan)
        label_dict = adata.var['highly_variable'].to_dict()
        adata2.var['highly_variable'] = adata2.var.index.map(label_dict).fillna(False)
        #adata2.write(adata_path + ".HVGs")
        print(f"Added HVGs")

        adata=adata2.copy()
        del adata2
        gc.collect()  # Force garbage collection

    else:
        adata=sc.read_h5ad(adata_path + ".HVGs")
try:
    adata.X = adata.layers["counts"].copy()
except:
    adata.layers["counts"] = adata.X.copy()
    print(adata.X[:6,:6].A)
    print("set adata.X to be above as no counts layer")
    # del(adata.layers["counts"])


best_HVG_BATCH_MINIMUM = None
closest_hvg_number = None
closest_difference = float('inf')
for HVG_BATCH_MINIMUM in [5,10,20,30, 50, 60, 70, 90, 100, 110, 120,130, 150,160,170,180, 190,200,210, 220, 240,250, 260, 270, 290, 300, 350, 375, 400, 450,470,480,490, 500,510,520,530, 550, 600, 700, 750, 800, 850, 900, 950, 1000, 1050, 1100, 1200]:
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
try:
    print(f"{hvg_number} selected -> {adata.shape}. batches {HVG_BATCH_MINIMUM}")
except:
    print(f"{hvg_number} selected -> {adata.shape}. batches error")


try:
    print(best_HVG_BATCH_MINIMUM)
except:
    print("no hvg minimum")

adata=adata.copy()
if run_scanvi==True:
    print("RUN SCANVI")
    def run_scvi(adata_hvg,  hvg_number , max_epochs, batch_size_vae, N_LATENT=10, N_LAYERS=1):
        DISPERSION =  'gene-batch'
        try:
            details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
        except:
            details="missing"
        adata_save_name = 'umap_' + details +"__1"
        print(adata_save_name)
        scvi.model.SCANVI.setup_anndata(adata_hvg, 
                                   #layer="counts",
                                   #categorical_covariate_keys=["tech"],
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
    if run_mrvi==False:
        print("RUN scvi")
        def run_scvi(adata_hvg,  hvg_number , max_epochs,  batch_size_vae, N_LATENT=10, N_LAYERS=1):
            DISPERSION = 'gene-batch'
            try:
                details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
            except:
                details="missingdetails"
            adata_save_name = 'umap_' + details +"__1"
            print(adata_save_name)
            scvi.model.SCVI.setup_anndata(adata_hvg, 
                                          #layer="counts",
                                            batch_key=HVG_BATCH_KEY,
                                            #                                labels_key="broad_annotation",
                                            # unlabeled_category="New/unlabelled/excluded"
                                           )
            model = scvi.model.SCVI(adata_hvg, 
                            dispersion=DISPERSION,
                            n_latent = N_LATENT, 
                            n_layers = N_LAYERS,
                           )
            model.train(max_epochs=max_epochs,             
                        early_stopping=True,
                        accelerator='gpu',
                       early_stopping_patience=5, #use_gpu =True, 
                       batch_size=batch_size_vae)
            print("model trained")
            return adata_hvg, model
    elif run_mrvi==True:
        print("RUN MRVI")
        sample_key = "info_id6"
        adata=adata[adata.obs["Site_status"]!="Postrx"]
        adata.obs["Site_status"] = adata.obs["Site_status"].cat.remove_unused_categories()
        def run_scvi(adata_hvg,  hvg_number , max_epochs,  batch_size_vae, N_LATENT=10, N_LAYERS=1):
            DISPERSION = 'gene'
            try:
                details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
            except:
                details="missingdetails"
            adata_save_name = 'umap_' + details +"__1"
            print(adata_save_name)
            MRVI.setup_anndata(adata_hvg, 
                                      # layer="counts",
                                           sample_key=sample_key,
                                            batch_key="dataset_id"
                                            #                                labels_key="broad_annotation",
                                            # unlabeled_category="New/unlabelled/excluded"
                                           )
            model = MRVI(adata_hvg, 
                            #dispersion=DISPERSION,
                            #n_latent = N_LATENT, 
                            #n_layers = N_LAYERS,
                           )
            model.train(max_epochs=max_epochs,             
                        early_stopping=True,
                       early_stopping_patience=5, #use_gpu =True, 
                        accelerator="gpu",
                       batch_size=batch_size_vae)
            print("model trained")
            return adata_hvg, model
adata, model_test = run_scvi(adata, 
                     hvg_number=hvg_number, 
                              max_epochs=MAX_EPOCHS, 
                          #    HVG_BATCH_MINIMUM=HVG_BATCH_MINIMUM, 
                          batch_size_vae=512,
                              #CATEGORICAL_COV=CAT_COVS,
                             # CONTINUOUS_COV=[],
                              #  HVG_BATCH_KEY=HVG_BATCH_KEY,
                                           N_LATENT=N_LATENT,
                                          N_LAYERS=N_LAYERS )

if run_scanvi==True:
    BASE_DIR = f'/nfs/team298/ls34/new_disease_atlas/model_XENIUMLABELT_{DATASET}_{hvg_number}/'

    model_test.save(BASE_DIR,
                    save_anndata=True,
                     overwrite=True)
else:
    BASE_DIR = f'/nfs/team298/ls34/new_disease_atlas/model_XENIUMLABELT_{DATASET}_{hvg_number}/'
    model_test.save(BASE_DIR,
                        save_anndata=True,
                     overwrite=True)




print(f"trained. now re-load adata: {adata_path}")
adata=sc.read_h5ad(adata_path)


print("mrvi false")
latent = model_test.get_latent_representation() 
adata.obsm["X_scvi"] = latent
neighbor_id = f"neighbor_{NEIGHBOR}" 

print("start neighbours")
sc.pp.neighbors(adata, use_rep = 'X_scvi', metric = "euclidean", n_neighbors=NEIGHBOR,key_added=neighbor_id)
print("neighbours done")

print("start umap")
sc.tl.umap(adata, min_dist=MIN_DIST, neighbors_key =neighbor_id ) 
print("finished umap")



# #adata.write(BASE_DIR + 'adata_all_scvi7.h5ad.countsonly')# compression="gzip")
# #print("saved", BASE_DIR + 'adata_all_scvi7.h5ad.countsonly')
# leidenres_list = [0.1]
# for leidenres in leidenres_list:
#     print("###", leidenres)
#     leiden_id = "leiden_res" + str(leidenres) # gayoso 1.2
#     sc.tl.leiden(adata, resolution=leidenres, key_added=leiden_id, neighbors_key=neighbor_id)
# print("prep save")
adata.write(BASE_DIR + 'adata_all_scvi7.h5ad.countsonly.clustered')# compression="gzip")
print(BASE_DIR + 'adata_all_scvi7.h5ad.countsonly.clustered')



# leidenres_list = [1]
# for leidenres in leidenres_list:
#     print("###", leidenres)
#     leiden_id = "leiden_res" + str(leidenres) # gayoso 1.2
#     sc.tl.leiden(adata, resolution=leidenres, key_added=leiden_id, neighbors_key=neighbor_id)
# print("prep save")
# #adata.write(BASE_DIR + 'adata_all_scvi5.h5ad', compression="gzip")
# #adata.X = adata.layers["counts"].copy()
# #del(adata.layers["counts"])
# adata.write(BASE_DIR + 'adata_all_scvi7.h5ad.countsonly.clustered_leidenres1')# compression="gzip")


# print("Saved to")
# print(BASE_DIR + 'adata_all_scvi7.h5ad.countsonly.clustered_leidenres1')
