#!/usr/bin/env python
# coding: utf-8
# %%

import os, sys
import random
import warnings
import logging
from datetime import datetime
# import gdown

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import seaborn as sns
import squidpy as sq
from matplotlib import gridspec
from sklearn.preprocessing import MinMaxScaler
from re import sub
import numpy as np
import pickle

from nichecompass.models import NicheCompass
from nichecompass.utils import (add_gps_from_gp_dict_to_adata,
                                create_new_color_dict,
                                compute_communication_gp_network,
                                visualize_communication_gp_network,
                                extract_gp_dict_from_mebocost_ms_interactions,
                                #extract_gp_dict_from_mebocost_es_interactions,
                                extract_gp_dict_from_nichenet_lrt_interactions,
                                extract_gp_dict_from_omnipath_lr_interactions,
                                #filter_and_combine_gp_dict_gps,
                                filter_and_combine_gp_dict_gps_v2,
                                generate_enriched_gp_info_plots)


"""
https://github.com/Lotfollahi-lab/nichecompass
here we use version 0.3.0
"""

"""
A full tutorial is available here:
https://nichecompass.readthedocs.io/en/latest/tutorials/notebooks/mouse_cns_spatial_reference_mapping.html

Here we illustrate use of our core dataset as a reference for mapping the validation data
"""
"""
make sure inside this path, you have the folders gene_annotations 
and gene_programs with the files
(available from https://github.com/Lotfollahi-lab/nichecompass/tree/main/data)
"""
handle='/lustre/scratch124/cellgen/haniffa/projects/developmental_fibroblasts/nobackup_output/nichecompasss/nichecompass/' 


# # Set up reference and query (important part)


"""
Load adata prepared in tutorial
"""
# no 2 = 1024 genes
# 2 = 2048 genes
ADATA_PATH="/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/nobackup_data/adata_core_plus_time_fornc_refquery2.h5ad.svg" #'/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_new.h5ad.final.filtered.svg'
#'/lustre/scratch126/cellgen/lotfollahi/ls34/nemo/adata_all.h5ad.clustered.clustered10.good.prenichecompass.svg'

#'/lustre/scratch124/cellgen/haniffa/projects/beacon/nobackup_data/adata_mintflownonresponse.h5ad.labelled.svg'
adata_vis=sc.read_h5ad(ADATA_PATH)  

# %%


# %%
n_epochs = 100
 
### Dataset ###
n_neighbors = 8 # for constructing spatial connectivities
n_sampled_neighbors = 4 # for model training
species = "human" # assume human as reference is human skin
spatial_key = "spatial"
mapping_entity_key = "mapping_entity"
dataset = f"XeniumTIME_REFQ"






"""
SPLIT INTO REFERENCE AND QUERY - completed in preparation tutorial
"""
 
query_batches = list(adata_vis[adata_vis.obs["batch_nc"]=="query"].obs["sample"].unique())
reference_batches = list(adata_vis[adata_vis.obs["batch_nc"]=="reference"].obs["sample"].unique())


dryrun=False
def read_and_qc(sample_name, wtsi, path='rawdata.h5ad'):
    """ This function reads anndata object.
    It also calculates QC metrics. Modify this function if required by your workflow.
    """
    print(path)
    adata = sc.read_h5ad(path)

    adata.uns['spatial'][sample_name] = adata.uns['spatial'].pop(list(adata.uns['spatial'])[0])
    adata.obs['label'] = list(adata.uns['spatial'])[0]
    adata.obs['WTSI_ID'] = wtsi

    # fix TypeError when read in obsm
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)
    # Calculate QC metrics
    from scipy.sparse import csr_matrix
    
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['MT'] = [gene.startswith('MT-') for gene in adata.var_names]
    # adata.obs['mt_frac'] = adata[:, adata.var['MT'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    
    # add sample name to obs names
    # adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    # adata.obs_names = adata.obs["sample"] \
    #                       + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    return adata


# %%


### Model ###
# AnnData keys
counts_key = "counts"
adj_key = "spatial_connectivities"
cat_covariates_keys = ["sample"]
gp_names_key = "nichecompass_gp_names"
active_gp_names_key = "nichecompass_active_gp_names"
gp_targets_mask_key = "nichecompass_gp_targets"
gp_targets_categories_mask_key = "nichecompass_gp_targets_categories"
gp_sources_mask_key = "nichecompass_gp_sources"
gp_sources_categories_mask_key = "nichecompass_gp_sources_categories"
latent_key = "nichecompass_latent"

# Architecture
cat_covariates_embeds_injection = ["gene_expr_decoder"]
cat_covariates_embeds_nums = [len(reference_batches) + len(query_batches)] ## number samples
cat_covariates_no_edges = [True]
conv_layer_encoder = "gatv2conv" # change to "gatv2conv" if enough compute and memory
active_gp_thresh_ratio = 0.01

# Trainer
n_epochs_all_gps = 25
lr = 0.001
lambda_edge_recon = 500000.
lambda_gene_expr_recon = 300.
lambda_l1_masked = 0. # increase if gene selection desired
lambda_l1_addon = 100.
edge_batch_size = 1024 # increase if more memory available
use_cuda_if_available = True
 

### Analysis ###
cell_type_key = "Annotation"
latent_leiden_resolution = 1
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
sample_key = "sample"
spot_size = 250
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"



warnings.filterwarnings("ignore")
# Get time of notebook execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%Y%m%d_%H%M%S")
current_timestamp += dataset  ## Change this for your own project label
current_timestamp


# %%


# Define paths
ga_data_folder_path = f"{handle}/data/gene_annotations"
gp_data_folder_path = f"{handle}/data/gene_programs"
so_data_folder_path = f"{handle}/data/spatial_omics"
omnipath_lr_network_file_path = f"{gp_data_folder_path}/omnipath_lr_network.csv"
collectri_tf_network_file_path = f"{gp_data_folder_path}/collectri_tf_network_{species}.csv"
nichenet_lr_network_file_path = f"{gp_data_folder_path}/nichenet_lr_network_v2_{species}.csv"
nichenet_ligand_target_matrix_file_path = f"{gp_data_folder_path}/nichenet_ligand_target_matrix_v2_{species}.csv"
mebocost_enzyme_sensor_interactions_folder_path = f"{gp_data_folder_path}/metabolite_enzyme_sensor_gps"
gene_orthologs_mapping_file_path = f"{ga_data_folder_path}/human_mouse_gene_orthologs.csv"
artifacts_folder_path = f"{handle}/artifacts"
model_folder_path = f"{artifacts_folder_path}/spatial_reference_mapping/{current_timestamp}/model"
figure_folder_path = f"{artifacts_folder_path}/spatial_reference_mapping/{current_timestamp}/figures"


# %%


if dryrun != True:
    os.makedirs(ga_data_folder_path, exist_ok=True)
    os.makedirs(gp_data_folder_path, exist_ok=True)
    os.makedirs(model_folder_path, exist_ok=True)
    os.makedirs(figure_folder_path, exist_ok=True)
    os.makedirs(so_data_folder_path, exist_ok=True)



### This part connects and retrieve information from database, so re-run if unfortunately the server is too busy and you get error


import omnipath as op

def extract_gp_dict_from_omnipath_lr_interactions(
        species =  "human",
        min_curation_effort: int=2,
        load_from_disk: bool=False,
        save_to_disk: bool=False,
        lr_network_file_path: str="../data/gene_programs/" \
                                            "omnipath_lr_network.csv",
        gene_orthologs_mapping_file_path: str="../data/gene_" \
                                                        "annotations/human_" \
                                                        "mouse_gene_orthologs.csv",
        plot_gp_gene_count_distributions: bool=True,
        gp_gene_count_distributions_save_path: str=None) -> dict:
    """
    Retrieve 724 human ligand-receptor interactions from OmniPath and extract
    them into a gene program dictionary. OmniPath is a database of molecular
    biology prior knowledge that combines intercellular communication data from
    many different resources (all resources for intercellular communication
    included in OmniPath can be queried via
    ´op.requests.Intercell.resources()´). If ´species´ is ´mouse´, orthologs
    from human interactions are returned.

    Parts of the implementation are inspired by 
    https://workflows.omnipathdb.org/intercell-networks-py.html (01.10.2022).

    Parameters
    ----------
    species:
        Species for which the gene programs will be extracted. The default is
        human. Human genes are mapped to mouse orthologs using a mapping file.
        NicheCompass contains a default mapping file stored under
        "<root>/data/gene_annotations/human_mouse_gene_orthologs.csv", which was
        created with Ensembl BioMart
        (http://www.ensembl.org/info/data/biomart/index.html).
    min_curation_effort: 
        Indicates how many times an interaction has to be described in a 
        paper and mentioned in a database to be included in the retrieval.
    load_from_disk:
        If ´True´, the OmniPath ligand receptor interactions will be loaded from
        disk instead of from the OmniPath library.
    save_to_disk:
        If ´True´, the OmniPath ligand receptor interactions will additionally 
        be stored on disk. Only applies if ´load_from_disk´ is ´False´.
    lr_network_file_path:
        Path of the file where the OmniPath ligand receptor interactions will be
        stored (if ´save_to_disk´ is ´True´) or loaded from (if ´load_from_disk´
        is ´True´).
    gene_orthologs_mapping_file_path:
        Path of the file where the gene orthologs mapping is stored if species
        is ´mouse´.
    plot_gp_gene_count_distributions:
        If ´True´, display the distribution of gene programs per number of
        source and target genes.
    gp_gene_count_distributions_save_path:
        Path of the file where the gene program gene count distribution plot
        will be saved if ´plot_gp_gene_count_distributions´ is ´True´.

    Returns
    ----------
    gp_dict:
        Nested dictionary containing the OmniPath ligand-receptor interaction
        gene programs with keys being gene program names and values being
        dictionaries with keys ´sources´, ´targets´, ´sources_categories´, and
        ´targets_categories´, where ´sources´ contains the OmniPath ligands,
        ´targets´ contains the OmniPath receptors, ´sources_categories´ contains
        the categories of the sources, and ´targets_categories´ contains
        the categories of the targets.
    """
    if not load_from_disk:
        # Define intercell_network categories to be retrieved (see
        # https://workflows.omnipathdb.org/intercell-networks-py.html,
        # https://omnipath.readthedocs.io/en/latest/api/omnipath.interactions.import_intercell_network.html#omnipath.interactions.import_intercell_network)
        intercell_df = op.interactions.import_intercell_network(
            include=["omnipath", "pathwayextra", "ligrecextra"])
        lr_interaction_df = intercell_df[
            (intercell_df["category_intercell_source"] == "ligand")
            & (intercell_df["category_intercell_target"] == "receptor")]
        if save_to_disk:
            lr_interaction_df.to_csv(lr_network_file_path, index=False)
    else:
        lr_interaction_df = pd.read_csv(lr_network_file_path, index_col=0)

    # Only keep curated interactions (see
    # https://r.omnipathdb.org/reference/filter_intercell_network.html)
    lr_interaction_df = lr_interaction_df[
        lr_interaction_df["curation_effort"] >= min_curation_effort]

    # Group receptors by ligands
    grouped_lr_interaction_df = lr_interaction_df.groupby(
        "genesymbol_intercell_source")["genesymbol_intercell_target"].agg(
            list).reset_index()

    # Resolve protein complexes into individual genes
    def compute_elementwise_func(lst, func):
        return [func(item) for item in lst]

    def resolve_protein_complexes(x):
        if not x:  # catches None, NaN, empty string
            return []
        if isinstance(x, float):  # just in case it's a NaN float
            return []
        if "COMPLEX:" not in x:
            return [x]
        else:
            # Example: split out complexes if your format is like "COMPLEX:A_B"
            return x.replace("COMPLEX:", "").split("_")

    grouped_lr_interaction_df["sources"] = grouped_lr_interaction_df[
        "genesymbol_intercell_source"].apply(
            lambda x: list(set(resolve_protein_complexes(x))))
    
    
    
    grouped_lr_interaction_df["sources_categories"] = grouped_lr_interaction_df[
        "sources"].apply(lambda x: ["ligand"] * len(x))
    grouped_lr_interaction_df["targets"] = grouped_lr_interaction_df[
        "genesymbol_intercell_target"].apply(
            lambda x: list(set([element for sublist in compute_elementwise_func(x, resolve_protein_complexes) for element in sublist])))
    grouped_lr_interaction_df["targets_categories"] = grouped_lr_interaction_df[
        "targets"].apply(lambda x: ["receptor"] * len(x))
    

    #Extract gene programs and store in nested dict
    gp_dict = {}
    for _, row in grouped_lr_interaction_df.iterrows():
        gp_dict[row["genesymbol_intercell_source"] +
                "_ligand_receptor_GP"] = {
                    "sources": row["sources"],
                    "targets": row["targets"],
                    "sources_categories": row["sources_categories"],
                    "targets_categories": row["targets_categories"]}
        
    if species == "mouse":
        # Create mapping df to map from human genes to mouse orthologs
        mapping_df = pd.read_csv(gene_orthologs_mapping_file_path)
        grouped_mapping_df = mapping_df.groupby(
            "Gene name")["Mouse gene name"].agg(list).reset_index()
        
        # Map all genes in the gene programs to their orthologs from the mapping
        # df or capitalize them if no orthologs are found (one human gene can
        # have multiple mouse orthologs)
        for _, gp in gp_dict.items():
            gp["sources"] = [element for nested_list_l1 in [
                list_element for nested_list_l2 in [
                    grouped_mapping_df[
                        grouped_mapping_df["Gene name"] == source][
                            "Mouse gene name"].values.tolist() if
                            source in grouped_mapping_df["Gene name"].values else
                            [[source.capitalize()]] for source in gp["sources"]]
                            for list_element in nested_list_l2]
                            for element in nested_list_l1]
            gp["targets"] = [element for nested_list_l1 in [
                list_element for nested_list_l2 in [
                    grouped_mapping_df[
                        grouped_mapping_df["Gene name"] == target][
                            "Mouse gene name"].values.tolist() if
                            target in grouped_mapping_df["Gene name"].values else
                            [[target.capitalize()]] for target in gp["targets"]]
                            for list_element in nested_list_l2]
                            for element in nested_list_l1]
            gp["sources_categories"] = ["ligand"] * len(gp["sources"])
            gp["targets_categories"] = ["receptor"] * len(gp["targets"])
    
    if plot_gp_gene_count_distributions:
        create_gp_gene_count_distribution_plots(
            gp_dict=gp_dict,
            gp_plot_label="OmniPath",
            save_path=gp_gene_count_distributions_save_path)
        
    return gp_dict
omnipath_gp_dict = extract_gp_dict_from_omnipath_lr_interactions(
    species="human",
    min_curation_effort=2,
    load_from_disk=False,
    save_to_disk=True,
    lr_network_file_path=omnipath_lr_network_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=False,
    gp_gene_count_distributions_save_path=f"/omnipath_gp_gene_count_distributions.svg")
#omnipath_gp_dict.head(5)

omnipath_gp_names = list(omnipath_gp_dict.keys())
random.shuffle(omnipath_gp_names)
omnipath_gp_name = omnipath_gp_names[0]
print(f"{omnipath_gp_name}: {omnipath_gp_dict[omnipath_gp_name]}")

nichenet_gp_dict = extract_gp_dict_from_nichenet_lrt_interactions(
    species=species,
    version="v2",
    keep_target_genes_ratio=1.,
    max_n_target_genes_per_gp=250,
    load_from_disk=False,
    save_to_disk=True,
    lr_network_file_path=nichenet_lr_network_file_path,
    ligand_target_matrix_file_path=nichenet_ligand_target_matrix_file_path,
    gene_orthologs_mapping_file_path=gene_orthologs_mapping_file_path,
    plot_gp_gene_count_distributions=True)

nichenet_gp_names = list(nichenet_gp_dict.keys())
random.shuffle(nichenet_gp_names)
nichenet_gp_name = nichenet_gp_names[0]
print(f"{nichenet_gp_name}: {nichenet_gp_dict[nichenet_gp_name]}")

mebocost_gp_dict = extract_gp_dict_from_mebocost_ms_interactions(
    dir_path=mebocost_enzyme_sensor_interactions_folder_path,
    species=species,
    plot_gp_gene_count_distributions=True)

mebocost_gp_names = list(mebocost_gp_dict.keys())
random.shuffle(mebocost_gp_names)
mebocost_gp_name = mebocost_gp_names[0]
print(f"{mebocost_gp_name}: {mebocost_gp_dict[mebocost_gp_name]}")
gp_dicts = [omnipath_gp_dict, nichenet_gp_dict, mebocost_gp_dict]
combined_gp_dict = filter_and_combine_gp_dict_gps_v2(
    gp_dicts,
    verbose=True)

print(f"Number of gene programs after filtering and combining: "
      f"{len(combined_gp_dict)}.")



### from jimmy lee 
def select_slide2(adata, s, s_col='sample'):
    """ This function selects the data for one slide from the spatial anndata object.
    :param adata: Anndata object with multiple spatial experiments
    :param s: name of selected experiment
    :param s_col: column in adata.obs listing experiment name for each location
    """
    slide = adata[adata.obs[s_col].isin([s]), :]
#     s_keys = list(slide.uns['spatial'].keys())
#     s_spatial = np.array(s_keys)[[s in k for k in s_keys]][0]
#     slide.uns['spatial'] = {s_spatial: slide.uns['spatial'][s_spatial]}
    return slide







# adata_reference=adata_vis.copy()
# #sc.read_h5ad(ADATA_PATH + ".reference")  



# # Add the GP dictionary as binary masks to the adata
# logging.info("add_gps_from_gp_dict_to_adata")
# add_gps_from_gp_dict_to_adata(
#     gp_dict=combined_gp_dict,
#     adata=adata_reference,
#     gp_targets_mask_key=gp_targets_mask_key,
#     gp_targets_categories_mask_key=gp_targets_categories_mask_key,
#     gp_sources_mask_key=gp_sources_mask_key,
#     gp_sources_categories_mask_key=gp_sources_categories_mask_key,
#     gp_names_key=gp_names_key,
#     min_genes_per_gp=2,
#     min_source_genes_per_gp=1,
#     min_target_genes_per_gp=1,
#     max_genes_per_gp=None,
#     max_source_genes_per_gp=None,
#     max_target_genes_per_gp=None)


# %%


def cast_adata_to_float32(adata, counts_key="counts"):
    adata.X = adata.X.astype("float32")                           # keeps sparse
    if counts_key in adata.layers:                                # fix the layer the model uses
        adata.layers[counts_key] = adata.layers[counts_key].astype("float32")
    for lyr in adata.layers:                                      # belt-and-braces
        adata.layers[lyr] = adata.layers[lyr].astype("float32")
    return adata
adata_reference = cast_adata_to_float32(adata_reference, counts_key=counts_key)




# %%


# # # Initialize model
# logging.info("NicheCompass")
# model = NicheCompass(adata_reference,
#                      counts_key=counts_key,
#                      adj_key=adj_key,
#                      cat_covariates_embeds_injection=cat_covariates_embeds_injection,
#                      cat_covariates_keys=cat_covariates_keys,
#                      cat_covariates_no_edges=cat_covariates_no_edges,
#                      cat_covariates_embeds_nums=cat_covariates_embeds_nums,
#                      gp_names_key=gp_names_key,
#                      active_gp_names_key=active_gp_names_key,
#                      gp_targets_mask_key=gp_targets_mask_key,
#                      gp_targets_categories_mask_key=gp_targets_categories_mask_key,
#                      gp_sources_mask_key=gp_sources_mask_key,
#                      gp_sources_categories_mask_key=gp_sources_categories_mask_key,
#                      latent_key=latent_key,
#                      conv_layer_encoder=conv_layer_encoder,
#                      active_gp_thresh_ratio=active_gp_thresh_ratio)



# # %%


# # Train model
# logging.info("model.train")
# model.train(n_epochs=n_epochs,
#             n_epochs_all_gps=n_epochs_all_gps,
#             lr=lr,
#             lambda_edge_recon=lambda_edge_recon,
#             lambda_gene_expr_recon=lambda_gene_expr_recon,
#             lambda_l1_masked=lambda_l1_masked,
#             edge_batch_size=edge_batch_size,
#             n_sampled_neighbors=n_sampled_neighbors,
#             use_cuda_if_available=use_cuda_if_available,
#             verbose=True)


# %%


# logging.info("sc.pp.neighbors")

# sc.pp.neighbors(model.adata,
#                 use_rep=latent_key,
#                 key_added=latent_key)

# # Compute UMAP embedding
# logging.info("sc.tl.umap")
# sc.tl.umap(model.adata,
#            neighbors_key=latent_key,
#           min_dist=0.1 
#           )




# # Save trained model
# logging.info("model.save")
# model.save(dir_path=f"{model_folder_path}/reference",
#            overwrite=True,
#            save_adata=True,
#            adata_file_name="adata.h5ad")

# logging.info("reference model saved!")


# %%


# adata_batch_list = []
# print("Processing query batches...")
# # for batch in query_batches:
# #     print(f"Processing batch {batch}...")
# #     print("Loading data...")
# #     adata_batch = sc.read_h5ad(
# #         f"{so_data_folder_path}/{dataset}_{batch}.h5ad")
# for batch in query_batches:
#     print(f"Processing batch {batch}...")
#     print("Loading data...")
#     adata_batch = select_slide2(adata_vis, batch)
#     print("Computing spatial neighborhood graph...\n")
#     # Compute (separate) spatial neighborhood graphs
#     sq.gr.spatial_neighbors(adata_batch,
#                             coord_type="generic",
#                             spatial_key=spatial_key,
#                             n_neighs=n_neighbors)
    
#     # Make adjacency matrix symmetric
#     adata_batch.obsp[adj_key] = (
#         adata_batch.obsp[adj_key].maximum(
#             adata_batch.obsp[adj_key].T))
#     adata_batch_list.append(adata_batch)
# adata_query = ad.concat(adata_batch_list, join="inner")
# try:
#     logging.info("query model pre-selection:", adata_query.shape)
# except:
#     1
# # adata_query.var["spatially_variable"] = adata_query.var_names.isin(sv_genes)
# # adata_query.var["keep_gene"] = adata_query.var["spatially_variable"]
# # adata_query = adata_query[:, adata_query.var["keep_gene"] == True]
# print(adata_query.shape)
# try:
#     logging.info("query model !", adata_query.shape)
# except:
#     1
# # %%


# for i in range(len(adata_batch_list)):
#     if i == 0: # first batch
#         after_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[0].shape[0],
#             (adata_query.shape[0] -
#             adata_batch_list[0].shape[0])))
#         batch_connectivities.append(sp.hstack(
#             (adata_batch_list[0].obsp[adj_key],
#             after_batch_connectivities_extension)))
#     elif i == (len(adata_batch_list) - 1): # last batch
#         before_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0],
#             (adata_query.shape[0] -
#             adata_batch_list[i].shape[0])))
#         batch_connectivities.append(sp.hstack(
#             (before_batch_connectivities_extension,
#             adata_batch_list[i].obsp[adj_key])))
#     else: # middle batches
#         before_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0], len_before_batch))
#         after_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0],
#             (adata_query.shape[0] -
#             adata_batch_list[i].shape[0] -
#             len_before_batch)))
#         batch_connectivities.append(sp.hstack(
#             (before_batch_connectivities_extension,
#             adata_batch_list[i].obsp[adj_key],
#             after_batch_connectivities_extension)))
#     len_before_batch += adata_batch_list[i].shape[0]
# adata_query.obsp[adj_key] = sp.vstack(batch_connectivities)


# # %%


# # Combine spatial neighborhood graphs as disconnected components
# batch_connectivities = []
# len_before_batch = 0

#adata_query.obs[mapping_entity_key] = "query"
#adata_query=sc.read_h5ad(ADATA_PATH + ".query")  

import scanpy as sc
ADATA_PATH = '/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/nobackup_data/adata_core_plus_RELAPSE_refquery2.h5ad.svg.query'
adata_query=sc.read_h5ad(ADATA_PATH)

# Add the GP dictionary as binary masks to the adata
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_gp_dict,
    adata=adata_query,
    gp_targets_mask_key=gp_targets_mask_key,
    gp_targets_categories_mask_key=gp_targets_categories_mask_key,
    gp_sources_mask_key=gp_sources_mask_key,
    gp_sources_categories_mask_key=gp_sources_categories_mask_key,
    gp_names_key=gp_names_key,
    min_genes_per_gp=2,
    min_source_genes_per_gp=1,
    min_target_genes_per_gp=1,
    max_genes_per_gp=None,
    max_source_genes_per_gp=None,
    max_target_genes_per_gp=None)

logging.info("query model saved!")


# # ref query

# %%


#load_timestamp = current_timestamp # uncomment if you trained the model in this notebook
#model_folder_path = f"{artifacts_folder_path}/spatial_reference_mapping/{load_timestamp}/model"

adata_query = cast_adata_to_float32(adata_query, counts_key=counts_key)

# Loading reference model with query data
print("Retrieving reference model...")

PATH = '/nfs/team298/ls34/adult_skin/final_models/nichecompass_adultskin_reference/'

model = NicheCompass.load(
    dir_path=f"{PATH}",
    adata=adata_query,
    adata_file_name="adata.h5ad",
    gp_names_key=gp_names_key,
    unfreeze_all_weights=False,
    unfreeze_cat_covariates_embedder_weights=True)



# Train model
model.train(n_epochs=n_epochs,
            n_epochs_all_gps=n_epochs_all_gps,
            lr=lr,
            lambda_edge_recon=lambda_edge_recon,
            lambda_gene_expr_recon=lambda_gene_expr_recon,
            lambda_l1_masked=lambda_l1_masked,
            edge_batch_size=edge_batch_size,
            n_sampled_neighbors=n_sampled_neighbors,
            use_cuda_if_available=use_cuda_if_available)


# Save trained model
model.save(dir_path=f"{model_folder_path}/query",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata.h5ad")

logging.info("final model saved!")


# %%


# # Integrate reference and query data
# adata_batch_list = [adata_reference, adata_query]
# adata_reference_query = ad.concat(adata_batch_list, join="inner")
# # adata_reference_query.var["spatially_variable"] = adata_reference_query.var_names.isin(sv_genes)
# # adata_reference_query.var["keep_gene"] = adata_reference_query.var["spatially_variable"]
# # adata_reference_query = adata_reference_query[:, adata_reference_query.var["keep_gene"] == True]
# logging.info("final model saved! " )
# print(adata_reference_query.shape)

# # Combine spatial neighborhood graphs as disconnected components
# batch_connectivities = []
# len_before_batch = 0
# for i in range(len(adata_batch_list)):
#     if i == 0: # first batch
#         after_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[0].shape[0],
#             (adata_reference_query.shape[0] -
#             adata_batch_list[0].shape[0])))
#         batch_connectivities.append(sp.hstack(
#             (adata_batch_list[0].obsp[adj_key],
#             after_batch_connectivities_extension)))
#     elif i == (len(adata_batch_list) - 1): # last batch
#         before_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0],
#             (adata_reference_query.shape[0] -
#             adata_batch_list[i].shape[0])))
#         batch_connectivities.append(sp.hstack(
#             (before_batch_connectivities_extension,
#             adata_batch_list[i].obsp[adj_key])))
#     else: # middle batches
#         before_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0], len_before_batch))
#         after_batch_connectivities_extension = sp.csr_matrix(
#             (adata_batch_list[i].shape[0],
#             (adata_reference_query.shape[0] -
#             adata_batch_list[i].shape[0] -
#             len_before_batch)))
#         batch_connectivities.append(sp.hstack(
#             (before_batch_connectivities_extension,
#             adata_batch_list[i].obsp[adj_key],
#             after_batch_connectivities_extension)))
#     len_before_batch += adata_batch_list[i].shape[0]
# adata_reference_query.obsp[adj_key] = sp.vstack(batch_connectivities)



# logging.info("integrated!")


# # %%


# # Add the GP dictionary as binary masks to the adata
# add_gps_from_gp_dict_to_adata(
#     gp_dict=combined_gp_dict,
#     adata=adata_reference_query,
#     gp_targets_mask_key=gp_targets_mask_key,
#     gp_targets_categories_mask_key=gp_targets_categories_mask_key,
#     gp_sources_mask_key=gp_sources_mask_key,
#     gp_sources_categories_mask_key=gp_sources_categories_mask_key,
#     gp_names_key=gp_names_key,
#     min_genes_per_gp=2,
#     min_source_genes_per_gp=1,
#     min_target_genes_per_gp=1,
#     max_genes_per_gp=None,
#     max_source_genes_per_gp=None,
#     max_target_genes_per_gp=None)


# # Load query model with the integrated data
# print("Retrieving query model...")
# model = NicheCompass.load(
#     dir_path=f"{model_folder_path}/query",
#     adata=adata_reference_query,
#     adata_file_name="adata.h5ad",
#     gp_names_key=gp_names_key)


# print("Computing reference query latent GP space...")
# model.adata.obsm[latent_key], _ = model.get_latent_representation(
#    adata=model.adata,
#    counts_key=counts_key,
#    adj_key=adj_key,
#    cat_covariates_keys=cat_covariates_keys,
#    only_active_gps=True,
#    return_mu_std=True,
#    node_batch_size=model.node_batch_size_)

# print("Computing active GPs...")
# model.adata.uns[model.active_gp_names_key_] = model.get_active_gps()

# # Compute latent neighbor graph
# sc.pp.neighbors(model.adata,
#                 use_rep=latent_key,
#                 key_added=latent_key)

# # Compute UMAP embedding
# sc.tl.umap(model.adata,
#            neighbors_key=latent_key)


# # Save model
# model.save(dir_path=f"{model_folder_path}/reference_query",
#            overwrite=True,
#            save_adata=True,
#            adata_file_name="adata.h5ad")
# print("final model saved!")
 
  

