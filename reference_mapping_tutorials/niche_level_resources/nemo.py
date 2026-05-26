import sys
sys.path.append("../../utils")

import logging
import os
import warnings
from typing import List, Tuple
import torch
from sklearn.metrics import silhouette_score
import gseapy as gp
import anndata as ad
import pandas as pd
import scanpy as sc
from datasets import load_from_disk
from matplotlib.colors import Normalize
import squidpy as sq
from tqdm import tqdm
import pickle
import gc


# from nichejepa.utils.evaluation import (get_top_gene_score,
#                                         get_top_gene_pairs
#                                         )
#from sklearn.metrics.cluster import adjusted_rand_score, normalized_mutual_info_score

import matplotlib.pyplot as plt

## This function takes a list of cell_ids and returns the corresponding dataset object.
def subset_by_cell_ids(dataset, cell_id_list):
    cell_id_set = set(cell_id_list)
    cell_ids = dataset['cell_id']
    indices = [i for i, cid in tqdm(enumerate(cell_ids), total=len(cell_ids), desc="Finding matching indices") if cid in cell_id_set]
    return dataset.select(indices)

"""
note that app should be installed from the spatial FM - see https://github.com/Lotfollahi-lab/
"""
from app.infer import (embed_dataset,
                       harmonize_adata,
                       tokenize_adata,
                       perturb_dataset,
                       harmonize_tokenize_embed_pipeline,
                       get_gene_embed,
                       get_average_gene_embed,
                       get_spatial_score,
                       get_emd_distance )

plt.rcParams['font.size'] = 5
plt.rcParams['text.usetex'] = False
plt.rcParams['svg.fonttype'] = 'none'

sc.set_figure_params(
    dpi=50,
    dpi_save=300,
    figsize=(3, 2),
    facecolor='white',
    fontsize=7
)
spot_size = 20

# Improve reproducibility of UMAP and Leiden
os.environ['NUMBA_CPU_NAME'] = 'generic'

# Where model is located
model_folder_path = '/nfs/team361/sb75/nichejepa-reproducibility/artifacts/models/18062025_082526_412'

# SET PATH TO SAVE TOKENIZED ADATAS
PATH_TO_TOKENIZED_ADATA = f'/lustre/scratch126/cellgen/lotfollahi/ls34/nemo/tokenized_adata_revision/'
# adata_path='/nfs/team298/ls34/adult_skin/final_adatas/adata_xenium_freeze_plus3d.h5ad.september.plusnewdata'
# adata_all=sc.read_h5ad(adata_path)
PATH = "/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_new.h5ad.final.filtered.svg"
adata_all=sc.read_h5ad(PATH)
#adata_all=adata_all[adata_all.obs["batch_nc"]=="query"]
# Attach ensembl IDs (gene names should be in .var_names)
file_path = '/nfs/team298/ls34/ensmbl_gene_5k.pkl'     # or the full path if you moved it

with open(file_path, "rb") as f:
    gene2ens = pickle.load(f)

if "ensembl_id" not in adata_all.var.columns:
    adata_all.var["ensembl_id"] = adata_all.var.index.map(gene2ens)
    adata_all.var["gene_name"] = adata_all.var.index
    adata_all.var.index=adata_all.var["gene_name"] 
    del(adata_all.var["gene_name"] )
adata_all.var.head()
# set cell type key
cell_type_key = 'lvl5_annotation'
emb_layer = None
batch_size = 128
pin_memory = False
num_workers = 12
agg_excluded_tokens = None
top_k = None
LENGTH=adata_all.obs["info_id6"].unique()
for i,SECTION in enumerate(adata_all.obs["info_id6"].unique()):
    print(SECTION)
    print(i, "/", len(LENGTH))
    dataset_name=SECTION
    save_dataset_path = PATH_TO_TOKENIZED_ADATA+ f'adata_{dataset_name}.h5ad'
    print(save_dataset_path)
    if not os.path.exists(save_dataset_path):
        adata=adata_all[adata_all.obs["info_id6"]==SECTION].copy()
        adata = harmonize_adata(adata)
        dataset = tokenize_adata(adata,
                             model_folder_path,
                             nproc = 4,
                             processing_mode = 'parallel')
        num_shards = 32
        dataset.save_to_disk(
                    save_dataset_path,
                    num_shards=num_shards)
    else:
        print("load")
        dataset = load_from_disk(save_dataset_path)

    output_embed = embed_dataset(
        dataset=dataset,
        model_folder_path=model_folder_path,
        emb_layer=emb_layer,
        agg_excluded_tokens=agg_excluded_tokens,
        top_k=top_k,
        batch_size=batch_size,
        pin_memory=pin_memory,
        num_workers=num_workers)
    # cell embedding (not used here)
    adata.obsm['cell_emb'] = output_embed['cell_emb']
    # niche embedding)
    adata.obsm['FM_niche_embedding'] = output_embed['neighborhood_emb']
    adata.write(save_dataset_path)
    del(adata)
    gc.collect()
 


ADATAS = []

for fname in sorted(os.listdir(PATH_TO_TOKENIZED_ADATA)):
    if not fname.endswith(".h5ad"):
        continue

    path = os.path.join(PATH_TO_TOKENIZED_ADATA, fname)

    ad = sc.read_h5ad(path)

    # Keep provenance (ABSOLUTELY do this)
    ad.obs["source_file"] = fname

    ADATAS.append(ad)

# Concatenate into a single AnnData
adata = sc.concat(
    ADATAS,
    axis=0,                # cells
    join="outer",          # union of genes
    label="batch",         # adds adata.obs["batch"]
    keys=[a.obs["source_file"].iloc[0] for a in ADATAS],
    index_unique="-",      # avoid obs index collisions
)

# Free memory
del ADATAS
gc.collect()


    
    
    
sc.pp.neighbors(adata,
                n_neighbors=15,
                use_rep='FM_niche_embedding',
                key_added='neighborhood')
sc.tl.umap(adata,
           neighbors_key='neighborhood')


sc.pl.umap(adata,
           neighbors_key='neighborhood',
           color="info_id6")
# Set dataset params
emb_key = 'neighborhood'
latent_leiden_resolution = 0.5
latent_cluster_key = f'{emb_key}_emb_leiden_res{str(latent_leiden_resolution).replace(".", "_")}'

sc.tl.leiden(adata,
             neighbors_key=emb_key,
             key_added=latent_cluster_key,
             resolution=latent_leiden_resolution)

adata_all=sc.read_h5ad(PATH)
