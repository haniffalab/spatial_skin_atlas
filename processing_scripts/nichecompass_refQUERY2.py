#!/usr/bin/env python
# coding: utf-8
# %%

# %%
n_svg = 3000
n_epochs = 100
 
### Dataset ###
n_neighbors = 8


dataset = f"XeniumBEACON_{n_svg}svg_n{n_neighbors}_REFQ2_FIG8"
ADATA_PATH= '/lustre/scratch126/cellgen/lotfollahi/ls34/nemo/adata_all.h5ad.clustered.clustered10.good.prenichecompass'
reference_batches =  ['Week 8 (resolved)_CE3-SKI-28-FO-4-S22-E1_a',
 'Lesional_CE5-SKI-28-FO-1-S22_replicate',
 'Week 8 (resolved)_CE3-SKI-28-FO-4-S22-E1_b',
 'BK39_Non-lesional Baseline',
 'Baseline_resolved_CE3-SKI-28-FO-1-S22-B1',
 'Week 8 (resolved)_CE5-SKI-27-FO-4-S22-E1',
 'BK23_Lesional Baseline',
 '3D_BK25_week12-D2',
 '3D_BK25_week12-D1orE1b',
 'Baseline_never_CE5-SKI-27-FO-2-S22_replicate',
 '3D_BK22_Lesional_baseline-C2',
 'BK51_Never Lesional',
 'BK20_Week 12',
 'BK30_Lesional Baseline',
 '3D_BK22_Lesional_baseline-A1',
 'Week 8 (resolved)_CE4-SKI-27-FO-3-S22-E2_b',
 'BK30_Week 12',
 '3D_BK22_Lesional_baseline-B1',
 'Lesional_Baseline_resolved_CE3-SKI-24-FO-1-S22_replicate',
 'Baseline_resolved_CE6-SKI-20-FO-1-S22-C2',
 'Baseline_never_CE5-SKI-27-FO-2-S22-C1',
 'BK49_wk8 Relapse',
 'Week 8 (resolved)_CE3-SKI-28-FO-4-S22_replicate',
 'BK21_Non-lesional Baseline',
 'BK18_Week 12',
 'BK24_Week 12',
 'Baseline_never_CE3-SKI-28-FO-2-S22_replicate',
 '3D_BK25_week12-B2',
 'BK39_Week 12',
 'BK22_Non-lesional Baseline',
 'BK27_Week 12',
 'Baseline_resolved_CE4-SKI-27-FO-1-S22-B2',
 'Baseline_resolved_CE6-SKI-20-FO-1-S22_replicate',
 'BK49_Past Lesional',
 'BK18_Non-lesional Baseline',
 '3D_BK22_Lesional_baseline-D1',
 'BK27_Lesional Baseline',
 'Baseline_resolved_CE6-SKI-28-FO-4-S22_replicate',
 'BK25_Lesional Baseline',
 'BK25_Week 12',
 'Lesional_CE6-SKI-28-FO-4-S22_replicate',
 'Baseline_resolved_CE4-SKI-27-FO-1-S22_replicate',
 'BK50_Never Lesional',
 'BK51_Past Lesional',
 'Baseline_resolved_CE5-SKI-27-FO-1-S22_replicate',
 'BK24_Non-lesional Baseline',
 'BK50_Past Lesional',
 'BK39_Lesional Baseline',
 'Lesional_CE5-SKI-28-FO-1-S22-A1',
 '3D_BK22_Lesional_baseline-B2',
 'BK25_Non-lesional Baseline',
 'BK20_Non-lesional Baseline',
 '3D_BK25_week12-D1orE1a',
 'Baseline_resolved_CE5-SKI-27-FO-1-S22-B1',
 '3D_BK25_week12-C1',
 'Week 8 (resolved)_CE6-SKI-28-FO-3-S22_replicate',
 'BK30_Day 14',
 '3D_BK22_Lesional_baseline-A2',
 'Week 8 (resolved)_CE6-SKI-28-FO-3-S22-E2',
 'BK18_Lesional Baseline',
 'Baseline_resolved_CE3-SKI-28-FO-1-S22_replicate',
 'Lesional_CE4-SKI-27-FO-4-S22-A2',
 'BK24_Lesional Baseline',
 'Baseline_never_CE4-SKI-21-FO-1-S22_replicate',
 'BK43_Never Lesional',
 'BK22_Lesional Baseline',
 '3D_BK25_week12-A2',
 'BK22_Week 12',
 'Week 8 (resolved)_CE4-SKI-27-FO-3-S22_replicate',
 'BK51_wk8 Relapse',
 'BK23_Non-lesional Baseline',
 'BK27_Non-lesional Baseline',
 'Week 8 (resolved)_CE4-SKI-27-FO-3-S22-E2_a',
 'Lesional_CE4-SKI-27-FO-4-S22_replicate',
 'Baseline_never_CE3-SKI-28-FO-2-S22-C1',
 'BK46_Never Lesional',
 'BK30_Non-lesional Baseline',
 'BK21_Week 12',
 'BK43_Past Lesional',
 'Week 8 (resolved)_CE5-SKI-27-FO-4-S22_replicate',
 'BK21_Lesional Baseline',
 'BK46_Past Lesional',
 '3D_BK25_week12-B1',
 'BK20_Lesional Baseline',
 'BK51_Past Lesional wk8 relaspe',
 'Baseline_never_CE4-SKI-21-FO-1-S22-C2',
 '3D_BK22_Lesional_baseline-D2',
 'Lesional_CE6-SKI-28-FO-4-S22-A1',
 'BK49_Past Lesional wk8 relaspe',
 'Baseline_resolved_CE6-SKI-28-FO-1-S22-B2',
 '3D_BK25_week12-C2',
 '3D_BK25_week12-A1',
 'Lesional_CE3-SKI-24-FO-1-S22-A1',
 'BK23_Week 12',
 'BK49_Never Lesional']
 
 

query_batches = ['BK21-SKI-27-FO-1-S8-A3',
 'BK51-SKI-27-FO-2-S9-B2',
 'BK23-SKI-27-FO-1-S8-B1',
 'BK21-SKI-27-FO1-S11-C1',
 'BK27-SKI-27-FO-5-S9-D1',
 'CE3-SKI-28-FO-1-S25-E1',
 'BK21-SKI-27-FO-1-S13-C2',
 'BK22-SKI-27-FO-2-S7-A1',
 'BK30-SKI-28-FO-1-S6-B2',
 'BK39-SKI-27-FO-1-S8-D2',
 'CE3-SKI-28-FO-1-S25-D1',
 'BK30-SKI-28-FO-1-S14-C2',
 'BK27-SKI-27-FO-1-S6-C1',
 'CE4-SKI-27-FO-1-S25-S29-S32',
 'BK51-SKI-27-FO-2-S4-S8-S6',
 'BK23-SKI-27-FO-5-S9-A2',
 'CE3-SKI-28-FO-1-S28-D2']

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


# %%


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



# ## Define Paramters



#'/nfs/team298/ls34/xenium_atlas/model_ALL_CLEAN_scanvi_ALL/adata_counts_integrated_final_colored.h5ad'
adata_path=ADATA_PATH

species = "human"
### FROM ADATA.OBS.SAMPLE

 

spatial_key = "spatial"
mapping_entity_key = "mapping_entity"

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
n_sampled_neighbors = 4
use_cuda_if_available = True
 

### Analysis ###
cell_type_key = "Annotation"
latent_leiden_resolution = 1
latent_cluster_key = f"latent_leiden_{str(latent_leiden_resolution)}"
sample_key = "sample"
spot_size = 250
differential_gp_test_results_key = "nichecompass_differential_gp_test_results"


# %%


warnings.filterwarnings("ignore")


# %%


# Get time of notebook execution for timestamping saved artifacts
now = datetime.now()
current_timestamp = now.strftime("%Y%m%d_%H%M%S")
current_timestamp += dataset  ## Change this for your own project label
current_timestamp


# %%




# %%

handle='/lustre/scratch124/cellgen/haniffa/projects/developmental_fibroblasts/nobackup_output/nichecompasss/nichecompass/' # make sure inside this path, you have the folders gene_annotations and gene_programs with the files


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


# ### This part connects and retrieve information from database, so re-run if unfortunately the server is too busy and you get error

# %%

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


























# ## Preparing Anndata

# %%


adata_vis=sc.read_h5ad(ADATA_PATH)  






# %%


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
for batch in reference_batches:
    print(f"Processing batch {batch}...")
    print("Loading data...")
    adata_batch = select_slide2(adata_vis, batch)


# %%


adata_batch_list = []
print("Processing reference batches...")
for batch in reference_batches:
    print(f"Processing batch {batch}...")
    print("Loading data...")
    adata_batch = select_slide2(adata_vis, batch)
    print(f"Size {adata_batch.shape}")
    print("Computing spatial neighborhood graph...\n")
    # Compute (separate) spatial neighborhood graphs
    logging.info("sq.gr.spatial_neighbors")
    #try:
    sq.gr.spatial_neighbors(adata_batch,
                                coord_type="generic",
                                spatial_key=spatial_key,
                                n_neighs=n_neighbors)
    #except:
    #    continue
    print(f"Spatial neighbours done ## {adata_batch.shape}")

    # Make adjacency matrix symmetric
    adata_batch.obsp[adj_key] = (
        adata_batch.obsp[adj_key].maximum(
            adata_batch.obsp[adj_key].T))
    adata_batch_list.append(adata_batch)
print("List made: ...")
for x in adata_batch_list:
    print(x.shape)
adata_reference = ad.concat(adata_batch_list, join="inner")




# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
            (adata_reference.shape[0] -
            adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[adj_key],
            after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_reference.shape[0] -
            adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_reference.shape[0] -
            adata_batch_list[i].shape[0] -
            len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key],
            after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
adata_reference.obsp[adj_key] = sp.vstack(batch_connectivities)
adata_reference.obs[mapping_entity_key] = "reference"


# ### Feature selection

# %%


# Filter spatially variable genes

logging.info("sq.gr.spatial_autocorr")
sq.gr.spatial_autocorr(adata_reference, mode="moran", genes=adata_vis.var_names)
sv_genes = adata_reference.uns["moranI"].index[:n_svg].tolist()
adata_reference.var["spatially_variable"] = adata_reference.var_names.isin(sv_genes)

adata_reference.var["keep_gene"] = adata_reference.var["spatially_variable"]
adata_reference = adata_reference[:, adata_reference.var["keep_gene"] == True]
print(f"Keeping {len(adata_reference.var_names)} spatially variable, highly "
      "variable or gene program relevant genes.")




# Add the GP dictionary as binary masks to the adata
logging.info("add_gps_from_gp_dict_to_adata")
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_gp_dict,
    adata=adata_reference,
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


# %%
def cast_adata_to_float32(adata, counts_key="counts"):
    adata.X = adata.X.astype("float32")                           # keeps sparse
    if counts_key in adata.layers:                                # fix the layer the model uses
        adata.layers[counts_key] = adata.layers[counts_key].astype("float32")
    for lyr in adata.layers:                                      # belt-and-braces
        adata.layers[lyr] = adata.layers[lyr].astype("float32")
    return adata
adata_reference = cast_adata_to_float32(adata_reference, counts_key=counts_key)



# # Initialize model
logging.info("NicheCompass")
model = NicheCompass(adata_reference,
                     counts_key=counts_key,
                     adj_key=adj_key,
                     cat_covariates_embeds_injection=cat_covariates_embeds_injection,
                     cat_covariates_keys=cat_covariates_keys,
                     cat_covariates_no_edges=cat_covariates_no_edges,
                     cat_covariates_embeds_nums=cat_covariates_embeds_nums,
                     gp_names_key=gp_names_key,
                     active_gp_names_key=active_gp_names_key,
                     gp_targets_mask_key=gp_targets_mask_key,
                     gp_targets_categories_mask_key=gp_targets_categories_mask_key,
                     gp_sources_mask_key=gp_sources_mask_key,
                     gp_sources_categories_mask_key=gp_sources_categories_mask_key,
                     latent_key=latent_key,
                     conv_layer_encoder=conv_layer_encoder,
                     active_gp_thresh_ratio=active_gp_thresh_ratio)

# Train model
logging.info("model.train")
model.train(n_epochs=n_epochs,
            n_epochs_all_gps=n_epochs_all_gps,
            lr=lr,
            lambda_edge_recon=lambda_edge_recon,
            lambda_gene_expr_recon=lambda_gene_expr_recon,
            lambda_l1_masked=lambda_l1_masked,
            edge_batch_size=edge_batch_size,
            n_sampled_neighbors=n_sampled_neighbors,
            use_cuda_if_available=use_cuda_if_available,
            verbose=True)



# %%
# Compute latent neighbor graph
logging.info("sc.pp.neighbors")

sc.pp.neighbors(model.adata,
                use_rep=latent_key,
                key_added=latent_key)

# Compute UMAP embedding
logging.info("sc.tl.umap")
sc.tl.umap(model.adata,
           neighbors_key=latent_key,
          min_dist=0.1 
          )


# %%
sc.tl.leiden(model.adata,
                resolution=latent_leiden_resolution,
                key_added=latent_cluster_key,
                neighbors_key=latent_key)



# Save trained model
logging.info("model.save")
model.save(dir_path=f"{model_folder_path}/reference",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata.h5ad")

logging.info("reference model saved!")



adata_batch_list = []
print("Processing query batches...")
# for batch in query_batches:
#     print(f"Processing batch {batch}...")
#     print("Loading data...")
#     adata_batch = sc.read_h5ad(
#         f"{so_data_folder_path}/{dataset}_{batch}.h5ad")
for batch in query_batches:
    print(f"Processing batch {batch}...")
    print("Loading data...")
    adata_batch = select_slide2(adata_vis, batch)
    print("Computing spatial neighborhood graph...\n")
    # Compute (separate) spatial neighborhood graphs
    sq.gr.spatial_neighbors(adata_batch,
                            coord_type="generic",
                            spatial_key=spatial_key,
                            n_neighs=n_neighbors)
    
    # Make adjacency matrix symmetric
    adata_batch.obsp[adj_key] = (
        adata_batch.obsp[adj_key].maximum(
            adata_batch.obsp[adj_key].T))
    adata_batch_list.append(adata_batch)
adata_query = ad.concat(adata_batch_list, join="inner")
adata_query.var["spatially_variable"] = adata_query.var_names.isin(sv_genes)
adata_query.var["keep_gene"] = adata_query.var["spatially_variable"]
adata_query = adata_query[:, adata_query.var["keep_gene"] == True]




# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
            (adata_query.shape[0] -
            adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[adj_key],
            after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_query.shape[0] -
            adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_query.shape[0] -
            adata_batch_list[i].shape[0] -
            len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key],
            after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
adata_query.obsp[adj_key] = sp.vstack(batch_connectivities)

adata_query.obs[mapping_entity_key] = "query"

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



load_timestamp = current_timestamp # uncomment if you trained the model in this notebook
model_folder_path = f"{artifacts_folder_path}/spatial_reference_mapping/{load_timestamp}/model"

adata_query = cast_adata_to_float32(adata_query, counts_key=counts_key)

# Loading reference model with query data
print("Retrieving reference model...")
model = NicheCompass.load(
    dir_path=f"{model_folder_path}/reference",
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

# Integrate reference and query data
adata_batch_list = [adata_reference, adata_query]
adata_reference_query = ad.concat(adata_batch_list, join="inner")
adata_reference_query.var["spatially_variable"] = adata_reference_query.var_names.isin(sv_genes)
adata_reference_query.var["keep_gene"] = adata_reference_query.var["spatially_variable"]
adata_reference_query = adata_reference_query[:, adata_reference_query.var["keep_gene"] == True]
logging.info("final model saved! " )
print(adata_reference_query.shape)

# Combine spatial neighborhood graphs as disconnected components
batch_connectivities = []
len_before_batch = 0
for i in range(len(adata_batch_list)):
    if i == 0: # first batch
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[0].shape[0],
            (adata_reference_query.shape[0] -
            adata_batch_list[0].shape[0])))
        batch_connectivities.append(sp.hstack(
            (adata_batch_list[0].obsp[adj_key],
            after_batch_connectivities_extension)))
    elif i == (len(adata_batch_list) - 1): # last batch
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_reference_query.shape[0] -
            adata_batch_list[i].shape[0])))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key])))
    else: # middle batches
        before_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0], len_before_batch))
        after_batch_connectivities_extension = sp.csr_matrix(
            (adata_batch_list[i].shape[0],
            (adata_reference_query.shape[0] -
            adata_batch_list[i].shape[0] -
            len_before_batch)))
        batch_connectivities.append(sp.hstack(
            (before_batch_connectivities_extension,
            adata_batch_list[i].obsp[adj_key],
            after_batch_connectivities_extension)))
    len_before_batch += adata_batch_list[i].shape[0]
adata_reference_query.obsp[adj_key] = sp.vstack(batch_connectivities)



logging.info("integrated!")

# Add the GP dictionary as binary masks to the adata
add_gps_from_gp_dict_to_adata(
    gp_dict=combined_gp_dict,
    adata=adata_reference_query,
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


# Load query model with the integrated data
print("Retrieving query model...")
model = NicheCompass.load(
    dir_path=f"{model_folder_path}/query",
    adata=adata_reference_query,
    adata_file_name="adata.h5ad",
    gp_names_key=gp_names_key)


print("Computing reference query latent GP space...")
model.adata.obsm[latent_key], _ = model.get_latent_representation(
   adata=model.adata,
   counts_key=counts_key,
   adj_key=adj_key,
   cat_covariates_keys=cat_covariates_keys,
   only_active_gps=True,
   return_mu_std=True,
   node_batch_size=model.node_batch_size_)

print("Computing active GPs...")
model.adata.uns[model.active_gp_names_key_] = model.get_active_gps()


# Compute latent neighbor graph
sc.pp.neighbors(model.adata,
                use_rep=latent_key,
                key_added=latent_key)

# Compute UMAP embedding
sc.tl.umap(model.adata,
           neighbors_key=latent_key)


# Save model
model.save(dir_path=f"{model_folder_path}/reference_query",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata.h5ad")
print("final model saved!")
 
# adata = sc.read_h5ad(ADATA_PATH)
# try:
#     adata.obsm[latent_key], _ = model.get_latent_representation()
#     sc.pp.neighbors(adata,
#                     use_rep=latent_key,
#                     key_added=latent_key)
#     sc.tl.umap(adata,
#                neighbors_key=latent_key,
#               min_dist=0.1
#               )
#     sc.tl.leiden(adata,
#                 resolution=latent_leiden_resolution,
#                 key_added=latent_cluster_key,
#                 neighbors_key=latent_key)
#     adata.write(f"{model_folder_path}/reference/adata_final_allgenes.h5ad")
# except:
#     try:
#         adata.uns[latent_key], _ = model.get_latent_representation()
#         sc.pp.neighbors(adata,
#                         use_rep=latent_key,
#                         key_added=latent_key)
#         sc.tl.umap(adata,
#                    neighbors_key=latent_key,
#                   min_dist=0.1
#                   )
#         sc.tl.leiden(adata,
#                     resolution=latent_leiden_resolution,
#                     key_added=latent_cluster_key,
#                     neighbors_key=latent_key)
#         adata.write(f"{model_folder_path}/reference/adata_final_allgenes.h5ad")
#     except:
#         print("failed save 2")





samples = model.adata.obs[sample_key].unique().tolist()


# %%


batch_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=cat_covariates_keys[0])


# Create plot of mapping entity annotations in physical and latent space
# groups = None
# save_fig = True
# file_path = f"{figure_folder_path}/" \
#             "batches_latent_physical_space.svg"

# fig = plt.figure(figsize=(12, 14))
# title = fig.suptitle(t=f"NicheCompass Batches " \
#                        "in Latent and Physical Space",
#                      y=0.96,
#                      x=0.55,
#                      fontsize=20)
# spec1 = gridspec.GridSpec(ncols=1,
#                           nrows=2,
#                           width_ratios=[1],
#                           height_ratios=[3, 2])
# spec2 = gridspec.GridSpec(ncols=len(samples),
#                           nrows=2,
#                           width_ratios=[1] * len(samples),
#                           height_ratios=[3, 2])
# axs = []
# axs.append(fig.add_subplot(spec1[0]))
# sc.pl.umap(adata=model.adata,
#            color=[mapping_entity_key],
#            groups=groups,
#            #palette=mapping_entity_colors,
#            title=f"Batches in Latent Space",
#            ax=axs[0],
#            show=False)
# for idx, sample in enumerate(samples):
#     axs.append(fig.add_subplot(spec2[len(samples) + idx]))
#     sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
#                   color=[mapping_entity_key],
#                   groups=groups,
#                   #palette=mapping_entity_colors,
#                   spot_size=spot_size,
#                   title=f"Batches in Physical Space \n"
#                         f"(Sample: {sample})",
#                   legend_loc=None,
#                   ax=axs[idx+1],
#                   show=False)

# # Create and position shared legend
# handles, labels = axs[0].get_legend_handles_labels()
# # lgd = fig.legend(handles,
# #                  labels,
# #                  loc="center left",
# #                  bbox_to_anchor=(0.98, 0.5))
# # axs[0].get_legend().remove()

# # Adjust, save and display plot
# plt.subplots_adjust(wspace=0.2, hspace=0.25)
# if save_fig:
#     fig.savefig(file_path,
#                 bbox_extra_artists=(lgd, title),
#                 bbox_inches="tight")
# plt.show()


# %%


batch_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=cat_covariates_keys[0])


# %%


# Create plot of batch annotations in physical and latent space
groups = None
save_fig = True
file_path = f"{figure_folder_path}/" \
            "batches_latent_physical_space.svg"

fig = plt.figure(figsize=(12, 14))
title = fig.suptitle(t=f"NicheCompass Batches " \
                       "in Latent and Physical Space",
                     y=0.96,
                     x=0.55,
                     fontsize=20)
spec1 = gridspec.GridSpec(ncols=1,
                          nrows=2,
                          width_ratios=[1],
                          height_ratios=[3, 2])
spec2 = gridspec.GridSpec(ncols=len(samples),
                          nrows=2,
                          width_ratios=[1] * len(samples),
                          height_ratios=[3, 2])
# axs = []
# axs.append(fig.add_subplot(spec1[0]))
# sc.pl.umap(adata=model.adata,
#            color=[cat_covariates_keys[0]],
#            groups=groups,
#            palette=batch_colors,
#            title=f"Batches in Latent Space",
#            ax=axs[0],
#            show=False)
# for idx, sample in enumerate(samples):
#     axs.append(fig.add_subplot(spec2[len(samples) + idx]))
#     sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
#                   color=[cat_covariates_keys[0]],
#                   groups=groups,
#                   palette=batch_colors,
#                   spot_size=spot_size,
#                   title=f"Batches in Physical Space \n"
#                         f"(Sample: {sample})",
#                   legend_loc=None,
#                   ax=axs[idx+1],
#                   show=False)

# # Create and position shared legend
# handles, labels = axs[0].get_legend_handles_labels()
# lgd = fig.legend(handles,
#                  labels,
#                  loc="center left",
#                  bbox_to_anchor=(0.98, 0.5))
# axs[0].get_legend().remove()

# # Adjust, save and display plot
# plt.subplots_adjust(wspace=0.2, hspace=0.25)
# if save_fig:
#     fig.savefig(file_path,
#                 bbox_extra_artists=(lgd, title),
#                 bbox_inches="tight")
# plt.show()


# %%


cell_type_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=cell_type_key)


# %%


# # Create plot of cell type annotations in physical and latent space
# groups = None
# save_fig = True
# file_path = f"{figure_folder_path}/" \
#             "cell_types_latent_physical_space.svg"

# fig = plt.figure(figsize=(12, 14))
# title = fig.suptitle(t=f"Cell Types " \
#                        "in Latent and Physical Space",
#                      y=0.96,
#                      x=0.55,
#                      fontsize=20)
# spec1 = gridspec.GridSpec(ncols=1,
#                           nrows=2,
#                           width_ratios=[1],
#                           height_ratios=[3, 2])
# spec2 = gridspec.GridSpec(ncols=len(samples),
#                           nrows=2,
#                           width_ratios=[1] * len(samples),
#                           height_ratios=[3, 2])
# axs = []
# axs.append(fig.add_subplot(spec1[0]))
# sc.pl.umap(adata=model.adata,
#            color=[cell_type_key],
#            groups=groups,palette=cell_type_colors,
#            title=f"Cell Types in Latent Space",
#            ax=axs[0],
#            show=False)
# for idx, sample in enumerate(samples):
#     axs.append(fig.add_subplot(spec2[len(samples) + idx]))
#     sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
#                   color=[cell_type_key],
#                   groups=groups,
#                   palette=cell_type_colors,
#                   spot_size=spot_size,
#                   title=f"Cell Types in Physical Space \n"
#                         f"(Sample: {sample})",
#                   legend_loc=None,
#                   ax=axs[idx+1],
#                   show=False)

# # Create and position shared legend
# handles, labels = axs[0].get_legend_handles_labels()
# lgd = fig.legend(handles,
#                  labels,
#                  loc="center left",
#                  bbox_to_anchor=(0.98, 0.5))
# axs[0].get_legend().remove()

# # Adjust, save and display plot
# plt.subplots_adjust(wspace=0.2, hspace=0.25)
# if save_fig:
#     fig.savefig(file_path,
#                 bbox_extra_artists=(lgd, title),
#                 bbox_inches="tight")
# plt.show()


# %%


# Compute latent Leiden clustering
sc.tl.leiden(adata=model.adata,
            resolution=latent_leiden_resolution,
            key_added=latent_cluster_key,
            neighbors_key=latent_key)


# %%


latent_cluster_colors = create_new_color_dict(
    adata=model.adata,
    cat_key=latent_cluster_key)


# %%


# # Create plot of latent cluster / niche annotations in physical and latent space
# groups = None # set this to a specific cluster for easy visualization, e.g. ["0"]
# save_fig = True
# file_path = f"{figure_folder_path}/" \
#             f"res_{latent_leiden_resolution}_" \
#             "niches_latent_physical_space.svg"

# fig = plt.figure(figsize=(12, 14))
# title = fig.suptitle(t=f"NicheCompass Niches " \
#                        "in Latent and Physical Space",
#                      y=0.96,
#                      x=0.55,
#                      fontsize=20)
# spec1 = gridspec.GridSpec(ncols=1,
#                           nrows=2,
#                           width_ratios=[1],
#                           height_ratios=[3, 2])
# spec2 = gridspec.GridSpec(ncols=len(samples),
#                           nrows=2,
#                           width_ratios=[1] * len(samples),
#                           height_ratios=[3, 2])
# axs = []
# axs.append(fig.add_subplot(spec1[0]))
# sc.pl.umap(adata=model.adata,
#            color=[latent_cluster_key],
#            groups=groups,
#            palette=latent_cluster_colors,
#            title=f"Niches in Latent Space",
#            ax=axs[0],
#            show=False)
# for idx, sample in enumerate(samples):
#     axs.append(fig.add_subplot(spec2[len(samples) + idx]))
#     sc.pl.spatial(adata=model.adata[model.adata.obs[sample_key] == sample],
#                   color=[latent_cluster_key],
#                   groups=groups,
#                   palette=latent_cluster_colors,
#                   spot_size=spot_size,
#                   title=f"Niches in Physical Space \n"
#                         f"(Sample: {sample})",
#                   legend_loc=None,
#                   ax=axs[idx+1],
#                   show=False)

# # Create and position shared legend
# handles, labels = axs[0].get_legend_handles_labels()
# lgd = fig.legend(handles,
#                  labels,
#                  loc="center left",
#                  bbox_to_anchor=(0.98, 0.5))
# axs[0].get_legend().remove()

# # Adjust, save and display plot
# plt.subplots_adjust(wspace=0.2, hspace=0.25)
# if save_fig:
#     fig.savefig(file_path,
#                 bbox_extra_artists=(lgd, title),
#                 bbox_inches="tight")
# plt.show()


# %%


save_fig = True
file_path = f"{figure_folder_path}/" \
            f"res_{latent_leiden_resolution}_" \
            f"niche_composition_batches.svg"

df_counts = (model.adata.obs.groupby([latent_cluster_key, cat_covariates_keys[0]])
             .size().unstack())
df_counts.plot(kind="bar", stacked=True, figsize=(10,10))
legend = plt.legend(bbox_to_anchor=(1, 1), loc="upper left", prop={'size': 10})
legend.set_title("Batch Annotations", prop={'size': 10})
plt.title("Batch Composition of Niches")
plt.xlabel("Niche")
plt.ylabel("Cell Counts")
if save_fig:
    plt.savefig(file_path,
                bbox_extra_artists=(legend,),
                bbox_inches="tight")


# %%


save_fig = True
file_path = f"{figure_folder_path}/" \
            f"res_{latent_leiden_resolution}_" \
            f"niche_composition_cell_types.svg"

df_counts = (model.adata.obs.groupby([latent_cluster_key, cell_type_key])
             .size().unstack())
df_counts.plot(kind="bar", stacked=True, figsize=(10,10))
legend = plt.legend(bbox_to_anchor=(1, 1), loc="upper left", prop={'size': 10})
legend.set_title("Cell Type Annotations", prop={'size': 10})
plt.title("Cell Type Composition of Niches")
plt.xlabel("Niche")
plt.ylabel("Cell Counts")
if save_fig:
    plt.savefig(file_path,
                bbox_extra_artists=(legend,),
                bbox_inches="tight")


# %%


# Check number of active GPs
active_gps = model.get_active_gps()
print(f"Number of total gene programs: {len(model.adata.uns[gp_names_key])}.")
print(f"Number of active gene programs: {len(active_gps)}.")


# %%


# Display example active GPs
gp_summary_df = model.get_gp_summary()
gp_summary_df[gp_summary_df["gp_active"] == True].head()


# %%


# Set parameters for differential gp testing
selected_cats = ["0"]
comparison_cats = "rest"
title = f"NicheCompass Strongly Enriched Niche GPs"
log_bayes_factor_thresh = 2.3
save_fig = True
file_path = f"{figure_folder_path}/" \
            f"/log_bayes_factor_{log_bayes_factor_thresh}" \
             "_niches_enriched_gps_heatmap.svg"


# %%


# Run differential gp testing
enriched_gps = model.run_differential_gp_tests(
    cat_key=latent_cluster_key,
    selected_cats=selected_cats,
    comparison_cats=comparison_cats,
    log_bayes_factor_thresh=log_bayes_factor_thresh)
enriched_number = len(enriched_gps)
print(f"# of enriched gps (for 0 vs rest): {enriched_number}")

if enriched_number < 2:
    enriched_gps = model.run_differential_gp_tests(
        cat_key=latent_cluster_key,
        selected_cats=selected_cats[0],
        comparison_cats=comparison_cats,
        log_bayes_factor_thresh=log_bayes_factor_thresh)
    enriched_number = len(enriched_gps)
    print(f"v2 - # of enriched gps (for 0 vs rest): {enriched_number}")
    
if enriched_number < 2:
    enriched_gps = model.run_differential_gp_tests(
        cat_key=latent_cluster_key,
        selected_cats=None,
        comparison_cats=comparison_cats,
        log_bayes_factor_thresh=log_bayes_factor_thresh)
    enriched_number = len(enriched_gps)
    print(f"v3 - # of enriched gps (For none vs rest_: {enriched_number}")


# %%


# Results are stored in a df in the adata object
model.adata.uns[differential_gp_test_results_key]
#print(f"Columns: ")


# %%





# %%


import pickle
with open(f'{model_folder_path}/nichecompass.pkl', 'wb') as file:
    pickle.dump({'gp_summary_df': gp_summary_df,
                 'enriched_gps': enriched_gps,
                # 'network_df': network_df,
                 #'visualize_communication_gp_network': visualize_communication_gp_network,
                 # 'enriched_gp_summary_df': enriched_gp_summary_df,
                },
                file)


# %%


# Visualize GP activities of enriched GPs across niches
df = model.adata.obs[[latent_cluster_key] + enriched_gps].groupby(latent_cluster_key).mean()


# %%


import pickle
with open(f'{model_folder_path}/nichecompass.pkl', 'wb') as file:
    pickle.dump({'gp_summary_df': gp_summary_df,
                 'enriched_gps': enriched_gps,
                 'df': df
                # 'network_df': network_df,
                 #'visualize_communication_gp_network': visualize_communication_gp_network,
                 # 'enriched_gp_summary_df': enriched_gp_summary_df,
                },
                file)


# %%


# Visualize GP activities of enriched GPs across niches
df = model.adata.obs[[latent_cluster_key] + enriched_gps].groupby(latent_cluster_key).mean()

scaler = MinMaxScaler()
normalized_columns = scaler.fit_transform(df)
normalized_df = pd.DataFrame(normalized_columns, columns=df.columns)
normalized_df.index = df.index

plt.figure(figsize=(16, 8))  # Set the figure size
ax = sns.heatmap(normalized_df,
            cmap='viridis',
            annot=False,
            linewidths=0)
plt.xticks(rotation=45,
           fontsize=8,
           ha="right"
          )
plt.xlabel("Gene Programs", fontsize=16)
plt.savefig(f"{figure_folder_path}/enriched_gps_heatmap.svg",
            bbox_inches="tight")


# %%


# # Store gene program summary of enriched gene programs
# save_file = True
# file_path = f"{figure_folder_path}/" \
#             f"/log_bayes_factor_{log_bayes_factor_thresh}_" \
#             "niche_enriched_gps_summary.csv"

# gp_summary_cols = ["gp_name",
#                    "n_source_genes",
#                    "n_non_zero_source_genes",
#                    "n_target_genes",
#                    "n_non_zero_target_genes",
#                    "gp_source_genes",
#                    "gp_target_genes",
#                    "gp_source_genes_importances",
#                    "gp_target_genes_importances"]

# enriched_gp_summary_df = gp_summary_df[gp_summary_df["gp_name"].isin(enriched_gps)]
# cat_dtype = pd.CategoricalDtype(categories=enriched_gps, ordered=True)
# enriched_gp_summary_df.loc[:, "gp_name"] = enriched_gp_summary_df["gp_name"].astype(cat_dtype)
# enriched_gp_summary_df = enriched_gp_summary_df.sort_values(by="gp_name")
# enriched_gp_summary_df = enriched_gp_summary_df[gp_summary_cols]

# if save_file:
#     enriched_gp_summary_df.to_csv(f"{file_path}")
# else:
#     display(enriched_gp_summary_df)


# # %%


# plot_label = f"log_bayes_factor_{log_bayes_factor_thresh}_cluster_{selected_cats[0]}_vs_rest"
# save_figs = True

# generate_enriched_gp_info_plots(
#     plot_label=plot_label,
#     model=model,
#     sample_key=sample_key,
#     differential_gp_test_results_key=differential_gp_test_results_key,
#     cat_key=latent_cluster_key,
#     cat_palette=latent_cluster_colors,
#     n_top_enriched_gp_start_idx=0,
#     n_top_enriched_gp_end_idx=10,
#     feature_spaces=samples, # ["latent"]
#     n_top_genes_per_gp=3,
#     save_figs=save_figs,
#     figure_folder_path=f"{figure_folder_path}/",
#     spot_size=spot_size)


# # %%


# try:
#     gp_name = "CCL19_ligand_receptor_target_gene_GP"
#     network_df = compute_communication_gp_network(
#         gp_list=[gp_name],
#         model=model,
#         group_key=latent_cluster_key,
#         n_neighbors=n_neighbors)

#     visualize_communication_gp_network(
#         adata=model.adata,
#         network_df=network_df,
#         figsize=(16, 7),
#         cat_colors=latent_cluster_colors,
#         edge_type_colors=["#1f77b4"], 
#         cat_key=latent_cluster_key,
#         save=True,
#         save_path=f"{figure_folder_path}/gp_network_{gp_name}.svg",
#         )
# except:
#     print("didn't work - CCL19")


# # %%


# try:
#     gp_name = "CCL21_ligand_receptor_target_gene_GP"
#     network_df = compute_communication_gp_network(
#         gp_list=[gp_name],
#         model=model,
#         group_key=latent_cluster_key,
#         n_neighbors=n_neighbors)

#     visualize_communication_gp_network(
#         adata=model.adata,
#         network_df=network_df,
#         figsize=(16, 7),
#         cat_colors=latent_cluster_colors,
#         edge_type_colors=["#1f77b4"], 
#         cat_key=latent_cluster_key,
#         save=True,
#         save_path=f"{figure_folder_path}/gp_network_{gp_name}.svg",
#         )
# except:
#     print("didn't work - CCL21 GP")


# # %%


# try:
#     gp_name = "IGFBP7_ligand_receptor_target_gene_GP"
#     network_df = compute_communication_gp_network(
#         gp_list=[gp_name],
#         model=model,
#         group_key=latent_cluster_key,
#         n_neighbors=n_neighbors)

#     visualize_communication_gp_network(
#         adata=model.adata,
#         network_df=network_df,
#         figsize=(16, 7),
#         cat_colors=latent_cluster_colors,
#         edge_type_colors=["#1f77b4"], 
#         cat_key=latent_cluster_key,
#         save=True,
#         save_path=f"{figure_folder_path}/gp_network_{gp_name}.svg",
#         )
# except:
#     print("didn't work 0 IGFBP7 GP")


# # %%


# import pickle
# with open(f'{model_folder_path}/nichecompass2.pkl', 'wb') as file:
#     pickle.dump({'gp_summary_df': gp_summary_df,
#                  'enriched_gps': enriched_gps,
#                  'df': df,
#                  #'network_df': network_df,
#                  #'visualize_communication_gp_network': visualize_communication_gp_network,
#                  # 'enriched_gp_summary_df': enriched_gp_summary_df,
#                 },
#                 file)


# # %%


# # Save model
logging.info("model.save")
model.save(dir_path=f"{model_folder_path}/reference_query",
           overwrite=True,
           save_adata=True,
           adata_file_name="adata_final.h5ad")
#adata_reduced_genes = sc.read_h5ad(f"{model_folder_path}/reference_query/adata_final.h5ad")




