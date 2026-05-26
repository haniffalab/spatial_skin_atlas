from scvi.external import MRVI
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

adata_tomap_path = ### add path to your adata here. make sure 
                        # 1.adata.obs["lvl5_annotation"] is set to ""New/unlabelled/excluded"" for all cells
                        # 2. adata.X is counts data
BASE = #"/nfs/team298/ls34/new_disease_atlas/" # base dir to save output

adata_reference_path ='/nfs/team298/ls34/new_disease_atlas/model_XENIUMLABELT_final_scrna_paper_1917/adata_all_scvi7.h5ad.countsonly.clustered_leidenres1'


DATASET="KoreanMapping"
adata_path = f'/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/nobackup_output/adata_{DATASET}_and_Atlas.h5ad'
N_LATENT=30
N_LAYERS=2
NEIGHBOR=20
MIN_DIST=0.1
HVG_BATCH_KEY = "sample_id"
HVG_NUMBER = 6000
HVG_BATCH_MINIMUM=80
MAX_EPOCHS=30
run_scanvi=True
SCANVI_LABELS_KEY="lvl5_annotation"
SCANVI_UNLABELLED="New/unlabelled/excluded"
run_mrvi=False
prep_HVGS_done = False


MERGE_DATASETS=True

#adata_path = '/nfs/team298/ls34/new_disease_atlas/adata_Wound_and_Atlas.h5ad'
if not prep_HVGS_done:
    if DATASET=="nonlesional":
        adata_path = '/nfs/team298/ls34/new_disease_atlas/adata_nonlesional_markedstar.h5ad'
    else:
        if MERGE_DATASETS:
            # this is clean inflamm atlas data
            #adata=sc.read_h5ad(base_dir+'adata_inflamm_scanvi6_HVG8024_EPOCHS10.h5ad.v2.fordeconv.noHS.reintegratedSCANVI.updated10backup')

            # with doublets etx.
            adata=sc.read_h5ad(adata_reference_path)


            adata_tomap=sc.read_h5ad(adata_tomap_path)

            #adata.obs.index = adata.obs.index + "_Atlas"
            #adata_tomap.obs.index = adata_tomap.obs.index + "_Wound"
            adata.var_names_make_unique()
            adata_tomap.var_names_make_unique()
            adata = ad.concat([adata, adata_tomap], label="Mapping_status2", keys=["Atlas", DATASET])
            adata.write(adata_path)
            print(f"appended {DATASET}")
        else:
            adata=sc.read_h5ad(adata_path)

    

    


    #adata_path = '/nfs/team298/ls34/new_disease_atlas/adata_all_scvi5.countsonly_notstar.h5ad'
    #'/nfs/team298/ls34/new_disease_atlas/adata_all_scvi5.h5ad.countsonly' - with Tstar

    # ADATA_PATH_FOR_LABELS = "/nfs/team298/ls34/disease_atlas/mrvi/adata_inflamm_scanvi6_HVG8024_EPOCHS10.h5ad.v2.fordeconv.noHS.reintegratedSCANVI.updated10"
    # adata_labelled = sc.read_h5ad(ADATA_PATH_FOR_LABELS)
    # label_dict = adata_labelled.obs["lvl5_annotation_new3"].to_dict()
    # import pickle

    # Now `labels_dict` contains the loaded dictionary
    # save_path = "/nfs/team298/ls34/labels_lvl5_annotation.pkl"
    # with open(save_path, 'wb') as file:
    #    pickle.dump(label_dict, file)
    # adata.obs['lvl5_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")

    # adata.obs['lvl3_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")

    # adata.obs['lvl5_annotation'] = adata.obs.index.map(label_dict).fillna("New/unlabelled/excluded")
    # adata.layers["counts"]=adata.X.copy()
    # adata.write(adata_path)
    # print(f"updated labels and saved to {adata_path}")



    # print("save adata")

    # adata_i = adata[(adata.obs["leiden_res0.2"]=="22")|
    #                (adata.obs["leiden_res0.2"]=="23")|
    #                 (adata.obs["leiden_res0.2"]=="2")
    #                ]
    # adata_i.write(base_dir + 'adata_scvi4_removedjunk_LC_MigDC_only.h5ad')
    # adata_i = adata[adata.obs["leiden_res0.2"]=="7"]
    # adata_i.write(base_dir + 'adata_scvi4_Schwannrelated.h5ad')

    try:
        adata.X = adata.layers["counts"].copy()
    except:
        adata.layers["counts"] = adata.X.copy()
        print(adata.X[:6,:6].A)
        print("set adata.X to be above as no counts layer")
        # del(adata.layers["counts"])

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
    label_dict = adata.var['highly_variable_nbatches'].to_dict()
    label_dict2 = adata.var['highly_variable'].to_dict()

    adata=sc.read_h5ad(adata_path)
    adata.var['highly_variable_nbatches'] = adata.var.index.map(label_dict).fillna(np.nan)
    adata.var['highly_variable'] = adata.var.index.map(label_dict2).fillna(False)

else:
    # need to remove this
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
for HVG_BATCH_MINIMUM in list(np.arange(50, 501, 10)):
    var_genes_batch = adata.var.highly_variable_nbatches > HVG_BATCH_MINIMUM
    var_select = adata.var.highly_variable_nbatches >= HVG_BATCH_MINIMUM
    var_genes = var_select.index[var_select]
    hvg_number = len(var_genes)
    
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
adata = adata[:, var_select].copy()
print(f"{hvg_number} selected -> {adata.shape}")


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
                                   layer="counts",
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
    if run_mrvi==False:
        print("RUN scvi")
        sample_key = "sample_id"
        def run_scvi(adata_hvg,  hvg_number , max_epochs,  batch_size_vae, N_LATENT=10, N_LAYERS=1):
            DISPERSION = 'gene-batch'
            try:
                details = "hvg" + str(hvg_number) +   '_'.join(CATEGORICAL_COV) + '_'.join(CONTINUOUS_COV) +  "_maxepochs" + str(max_epochs) + "_nlatent" + str(N_LATENT)+"nlayers" + str(N_LAYERS) + "_BATCHKEY_" + HVG_BATCH_KEY.replace("_", "").lower() 
            except:
                details="missingdetails"
            adata_save_name = 'umap_' + details +"__1"
            print(adata_save_name)
            scvi.model.SCVI.setup_anndata(adata_hvg, 
                                          layer="counts",
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
        sample_key = "sample_id"
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
    BASE_DIR = BASE + f'/model_scanvi5_{DATASET}_{hvg_number}/'

    model_test.save(BASE_DIR,
                    save_anndata=True,
                     overwrite=True)
else:
    BASE_DIR = BASE + f'/model_scvi5_{DATASET}_{hvg_number}/'
    model_test.save(BASE_DIR,
                        save_anndata=True,
                     overwrite=True)


print(f"trained. now re-load adata: {adata_path}")
adata=sc.read_h5ad(adata_path)

# adata.layers["counts"] = adata.X.copy()
# sc.pp.normalize_total(adata, target_sum=1e4)
# sc.pp.log1p(adata)

if not run_mrvi==True:
    print("mrvi false")
    latent = model_test.get_latent_representation() 
    adata.obsm["X_scvi"] = latent
    neighbor_id = f"neighbor_{NEIGHBOR}"   
    print("start neighbours")
    sc.pp.neighbors(adata, use_rep = 'X_scvi', metric = "euclidean", n_neighbors=NEIGHBOR,key_added=neighbor_id)
    print("neighbours done")

    # if run_scanvi==True:
    #      model_test.save(f'/nfs/team298/ls34/disease_atlas/mrvi/v2_all_SCANVIversion_{hvg_number}',
    #                     save_anndata=True) 
    # else:
    #     model_test.save(f'/nfs/team298/ls34/disease_atlas/mrvi/v2_all_SCVIversion_{hvg_number}',
    #                         save_anndata=True)
    print("start umap")
    sc.tl.umap(adata, min_dist=MIN_DIST, neighbors_key =neighbor_id ) 
    #adata.write(f'/nfs/team298/ls34/disease_atlas/mrvi/adata_mrvi_all_scvi4_v2.h5ad')
    print("finished umap")
    leidenres_list = [0.1]
    for leidenres in leidenres_list:
        print("###", leidenres)
        leiden_id = "leiden_res" + str(leidenres)  
        sc.tl.leiden(adata, resolution=leidenres, key_added=leiden_id, neighbors_key=neighbor_id)
    print("prep save")
    adata.write(BASE_DIR + 'adata_all_scvi5.h5ad', compression="gzip")
    print("Saved to")
    print(BASE_DIR + 'adata_all_scvi5.h5ad')


#     cluster_counts = (
#         adata.obs.groupby([leiden_id, "lvl5_annotation"])
#         .size()
#         .reset_index(name="count")
#     )
#     cluster_counts = (
#         cluster_counts.sort_values([leiden_id, "count"], ascending=[True, False])
#         .groupby(leiden_id)
#         .head(2)  # Keep top 2 most common per cluster
#     )
#     def select_best_annotation(group):
#         top_annotations = group["lvl5_annotation"].tolist()
#         if top_annotations[0] == SCANVI_UNLABELLED and len(top_annotations) > 1:
#             return top_annotations[1]  # Use second most common
#         return top_annotations[0]  # Otherwise, use the most common

#     provisional_types = (
#         cluster_counts.groupby(leiden_id).apply(select_best_annotation).to_dict()
#     )
#     adata.obs["provisional_celltypes"] = adata.obs[leiden_id].map(provisional_types)

    

#     adata.write(BASE_DIR + 'adata_all_scvi5.h5ad', compression="gzip")
#     print("Saved to")
#     print(BASE_DIR + 'adata_all_scvi5.h5ad')

                    


    from datetime import datetime
    #sc.pp.subsample(adata,0.2)
    #adata.write(f'/nfs/team298/ls34/disease_atlas/mrvi/adata_all_scvi5_healthy.h5ad.countsonly.subsampled', compression="gzip")
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"Saved! Time: {timestamp}")

else:
    adata=adata[adata.obs["Site_status"]!="Postrx"]
    adata.obs["Site_status"] = adata.obs["Site_status"].cat.remove_unused_categories()

    latent = model_test.get_latent_representation() 
    adata.obsm["X_mrvi_u"] = latent
    latent_z = model_test.get_latent_representation(give_z=True) 
    adata.obsm["X_mrvi_z"] = latent_z
    neighbor=30
    neighbor_id = "neighbor_30_U"   
    sc.pp.neighbors(adata, use_rep = 'X_mrvi_z', metric = "euclidean", n_neighbors=neighbor,key_added="neighbor_30_Z")
    sc.pp.neighbors(adata, use_rep = 'X_mrvi_u', metric = "euclidean", n_neighbors=neighbor,key_added=neighbor_id)
    print("neighbours done")
    sc.pp.neighbors(adata, use_rep = 'X_mrvi_z', metric = "euclidean", n_neighbors=neighbor,key_added="neighbor_30_Z")
    mindist=0.2
    sc.tl.umap(adata, min_dist=mindist, neighbors_key =neighbor_id ) 
    leidenres_list = [0.5]
    for leidenres in leidenres_list:
        print("###", leidenres)
        leiden_id = "leiden_res" + str(leidenres) # gayoso 1.2
        sc.tl.leiden(adata, resolution=leidenres, key_added=leiden_id, neighbors_key=neighbor_id)
    print("prep save")

    dists = model_test.get_local_sample_distances(
    keep_cell=False, groupby="lvl3_annotation", batch_size=32
    )
    print("dists below")
    print(dists)
    try:
        d1 = dists.loc[{"lvl3_annotation_name": "Th"}].initial_clustering
    except:
        print("fail with: dists.loc[{lvl3_annotation_name: Th}].initial_clustering") 
    try:
        print(model_test.sample_info.columns)
    except:
        print("cant print model_test.sample_info.columns")
    try:
        print(model_test.sample_info["Site_status"].value_counts())
    except:
        print("cant print model_test.sample_info[status].value_counts")
    try:
        sample_cov_keys = ["Site_status"]
        model_test.sample_info["Site_status"] = model_test.sample_info["Site_status"].cat.reorder_categories(
            ["Nonlesional", "Lesional"]
        )  
        de_res = model_test.differential_expression(
            sample_cov_keys=sample_cov_keys,
            store_lfc=True,
        )
        print("de_Res done")
    except:
        print("faill sample_imnfo")
    try:
        print(de_res)
    except:
        print("cant print deres")
    try:
        adata.obs["DE_eff_size"] = de_res.effect_size.sel(covariate="Status_Lesional").values
        print("DE_EFF_SIZE added")
    except:
        print("fail Covid_DE_eff_size")
    try:
        print("da_res below")
        print(model_test.differential_abundance)
    except:
        print("could not print model_test.differential_abundance")
    try:
        sample_cov_keys = ["Site_status"]
        da_res = model_test.differential_abundance(sample_cov_keys=sample_cov_keys)
        covid_log_probs = da_res.Status_log_probs.loc[{"Site_status": "Lesional"}]
        healthy_log_probs = da_res.Status_log_probs.loc[{"Site_status": "Nonlesional"}]
        covid_healthy_log_prob_ratio = covid_log_probs - healthy_log_probs
        adata.obs["DA_lfc"] = covid_healthy_log_prob_ratio.values
    except:
        print(f"fail. da_res is {da_rest}")

    adata.write(f'/nfs/team298/ls34/disease_atlas/mrvi/adata_mrvi_all_scvi4.h5ad')
    from datetime import datetime
    now = datetime.now()
    timestamp = now.strftime("%Y-%m-%d %H:%M:%S")
    print(f"Saved! Time: {timestamp} -> /nfs/team298/ls34/disease_atlas/mrvi/adata_mrvi_all.h5ad")    
    model_test.save(f'/nfs/team298/ls34/disease_atlas/mrvi/scvi4_all_MRVIversion_{hvg_number}_new',
                   save_anndata=True) 
    #sc.settings.figdir="/nfs/team298/ls34/disease_atlas/mrvi/"

