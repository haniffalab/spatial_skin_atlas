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
import pickle

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

# DATASET="v2skin"f
out_dir = "/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/tmp/"
sc.settings.set_figure_params(dpi_save=300, facecolor="white", frameon=False, figsize=(20,20))
sc.settings.figdir = out_dir  


SAVE_PATH = '/nfs/team298/ls34/adult_skin/final_adatas/adata_combined_integrated2.h5ad'
BASE = '/lustre/scratch124/cellgen/haniffa/projects/adult_skin_v1/nobackup_output/subclustered_adatas_v2/'
adata=sc.read_h5ad(SAVE_PATH)

for x in adata.obs["coarse"].unique():
   # if x not in ["LTo", "YS_Mesothelium"]:
    adata_i = adata[adata.obs["coarse"]==x]

    N_NEIGHBORS = 20
    MIN_DIST = 0.1

    sc.pp.neighbors(adata_i, n_neighbors=N_NEIGHBORS, use_rep = "X_scvi", key_added=f"n_{N_NEIGHBORS}")
    print(f"neighbors done", N_NEIGHBORS)
    sc.tl.umap(adata_i, min_dist=MIN_DIST, neighbors_key =f"n_{N_NEIGHBORS}") 
    #sc.tl.leiden(adata_i, resolution=3, key_added='leiden_res3', neighbors_key="neighbor_20")
    sc.tl.leiden(adata_i, resolution=2, key_added='leiden_res2', neighbors_key="neighbor_20")
    x=x.replace("/", "_")

    adata_i.write(BASE + f'adata_subclustered_{x}.h5ad')

    sc.tl.paga(adata_i, groups="combined_anno",   neighbors_key=f"n_{N_NEIGHBORS}")
    sc.pl.paga(adata_i,
               color="combined_anno",
               #neighbors_key=f"n_{N_NEIGHBORS}",
               add_pos=True,    # ← this computes adata_i.uns['paga']['pos']
               show=True)
    sc.tl.umap(adata_i, init_pos='paga',   min_dist=0.3,
             neighbors_key = f"n_{N_NEIGHBORS}"
              )
    adata_i.write(BASE + f'adata_subclustered_{x}.h5ad.PAGA')





