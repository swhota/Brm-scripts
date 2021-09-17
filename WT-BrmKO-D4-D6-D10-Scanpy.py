import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import gseapy
import collections
import seaborn as sns
import bbknn
import scipy.stats
from scanpy.tools._utils import get_init_pos_from_paga as init

sc.set_figure_params(dpi=100, color_map='viridis')
sc.settings.verbosity = 2
sc.logging.print_versions()

class preprocess():
    def __init__(self, adata):
        self.adata = adata
    
    def high_exp_genes(self):
        print('Genes that yield the highest fraction of counts in each single cells, across all cells')
        sc.pl.highest_expr_genes(self.adata, n_top=20)
    
    def basic_filtering(self):
        sc.pp.filter_cells(self.adata, min_genes=200)
        sc.pp.filter_genes(self.adata, min_cells=3)
    
    def gene_percentage(self, gene_start, col_name):
        gene_boolean_list = self.adata.var_names.str.startswith(gene_start)
        
        self.adata.var[gene_start] = gene_boolean_list
        self.adata.uns[gene_start] = self.adata.var[self.adata.var[gene_start].isin([True])].index.tolist()

        print('Number of genes that start with ' + gene_start + ' {}'.format(collections.Counter(self.adata.var_names.str.startswith(gene_start))[True]))
        # for each cell compute fraction of counts in genes vs. all genes
        # the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
        self.adata.obs[col_name] = np.sum(
        self.adata[:, gene_boolean_list].X, axis=1).A1 / np.sum(self.adata.X, axis=1).A1
        
    
    def qc_metrics(self):
        self.gene_percentage('mt-', 'percent_mito')
        self.gene_percentage('Rpl', 'percent_rpl')
        self.gene_percentage('Rps', 'percent_rps')
        self.gene_percentage('Hba-', 'percent_hba')
        self.gene_percentage('Hbb-', 'percent_hbb')
      
        # add the total counts per cell as observations-annotation to adata
        self.adata.obs['n_counts'] = self.adata.X.sum(axis=1).A1
        print('Computed quality metrics')
        sc.pl.violin(self.adata, ['n_genes', 'n_counts', 'percent_mito',\
                                  'percent_rpl', 'percent_rps',\
                                 'percent_hba', 'percent_hbb'], jitter=0.4, multi_panel=True)
        
        sns.relplot(x="n_counts", y="n_genes", hue="niceorder",
            sizes=(40, 400), alpha=.05, palette="muted",
            height=6, data=self.adata.obs)
        
        sns.relplot(x="n_counts", y="percent_mito", hue="niceorder",
            sizes=(40, 400), alpha=.05, palette="muted",
            height=6, data=self.adata.obs)
        
        sns.relplot(x="n_counts", y="percent_rpl", hue="niceorder",
            sizes=(40, 400), alpha=.05, palette="muted",
            height=6, data=self.adata.obs)
        
        sns.relplot(x="n_counts", y="percent_rps", hue="niceorder",
            sizes=(40, 400), alpha=.05, palette="muted",
            height=6, data=self.adata.obs)
        
        return self.adata
    
    def filter_qc_data(self):
        self.adata = self.adata[self.adata.obs['n_genes'] < 2500, :]
        self.adata = self.adata[self.adata.obs['percent_mito'] < 0.05, :]
        return self.adata
    
    def normalize(self):
        sc.pp.normalize_per_cell(self.adata, counts_per_cell_after=1e4)
        sc.pp.log1p(self.adata)
    
    def high_var_genes(self):
#         n_top_genes=2000, batch_key='day'
        sc.pp.highly_variable_genes(self.adata,  n_top_genes = 2000, batch_key='diff_day')
        print('Highly variable genes')
        sc.pl.highly_variable_genes(self.adata)
        return self.adata
    
    def regress_vars_mito_and_rp(self):
        sc.pp.regress_out(self.adata, ['n_counts', 'percent_mito', 'percent_rpl', 'percent_rps'])
    
    def regress_vars_mito(self):
        sc.pp.regress_out(self.adata, ['n_counts', 'percent_mito'])
        
    def scale(self):
        # try to remove max value
        # default is 10
        sc.pp.scale(self.adata, max_value=10)

    def pca(self):
        sc.tl.pca(self.adata, svd_solver='arpack')
    
    def neighborhood_graph(self):
        # 10 and 40 are default
        sc.pp.neighbors(self.adata, n_neighbors=10, n_pcs=40)
    
    def bbknn_graph(self):
        bbknn.bbknn(self.adata, batch_key='geno', neighbors_within_batch=3, trim=5)
    
    def umap(self):
        sc.tl.umap(self.adata)
        
    def louvain(self):
        sc.tl.louvain(self.adata, flavor='vtraag')
    
    def leiden(self):
        sc.tl.leiden(self.adata, resolution=0.7)
    
    def run_up_to_high_var(self):
        self.high_exp_genes()
        self.basic_filtering()
        self.qc_metrics()
        self.adata = self.filter_qc_data()
        self.normalize()
        self.adata.raw = self.adata
        self.adata = self.high_var_genes()
        return self.adata
    
    def run(self, rp_regress=True):
        self.high_exp_genes()
        self.basic_filtering()
        self.qc_metrics()
        self.adata = self.filter_qc_data()
        self.normalize()
        self.adata.raw = self.adata
        self.adata = self.high_var_genes()
        self.adata = self.adata[:, self.adata.var['highly_variable']]
        if rp_regress:
            self.regress_vars_mito_and_rp()
        else:
            self.regress_vars_mito()
        self.scale()
        self.pca()
#         self.bbknn_graph()
        self.neighborhood_graph()
        self.umap()
        self.louvain()
        self.leiden()
        sc.pl.umap(self.adata, color=['percent_mito', 'percent_rpl', 'percent_rps', 'niceorder', 'louvain', 'leiden'])
        return self.adata

class PagaTrajectory():
    def __init__(self, adata, group):
        self.adata = adata
        self.group = group
    
    def build_paga_graph(self):
        sc.tl.paga(self.adata, groups=self.group)
    
    def view_course_graph(self):
        for meta in ['geno', 'day', 'niceorder', 'louvain']:
            sc.pl.umap(self.adata, color=meta)
            sc.pl.paga(self.adata, color=meta)
    
    def reconstruct_with_umap_embedding(self):
        sc.tl.umap(self.adata, init_pos=init(self.adata))
        print('UMAP reconstructed with PAGA coords')
        sc.pl.umap(self.adata, color=['niceorder'])
    
    def paga_compare(self):
        sc.pl.paga_compare(self.adata, color=['niceorder'], 
                           threshold=0.03, title='', right_margin=0.2,
                           size=10, edge_width_scale=0.5,
                           legend_fontsize=12, fontsize=12,
                           frameon=False, edges=True, save=True)
    
    def run(self):
        self.build_paga_graph()
        self.view_course_graph()
        self.reconstruct_with_umap_embedding()
#         self.paga_compare()
        return self.adata

def split_adata(brm_adata):
    # Low bmp condition
    low_bmp_adata = brm_adata[brm_adata.obs['niceorder'].isin([genotype for
                                                            genotype in 
                                                            list(set(brm_adata.obs['niceorder'])) if genotype.split('_')[-1] == 'lowBMP4'])]
    # high bmp condition
    high_bmp_adata = brm_adata[brm_adata.obs['niceorder'].isin([genotype for
                                                            genotype in 
                                                            list(set(brm_adata.obs['niceorder'])) if genotype.split('_')[-1] == 'highBMP4'])]
    return low_bmp_adata, high_bmp_adata

def color_tuple(adata, category):
    color_df = pd.DataFrame([list(color_vec) for color_vec in zip(list(adata.obs[category].cat.categories), adata.uns[category + '_colors'])], columns=[category, 'hex'])
    color_df = color_df.set_index(category)
    color_df['new_hex'] = 0
    return color_df

# Read in full dataset
brm_adata = sc.read('../Data/BrmKO_deep.h5ad')
low_bmp_adata, high_bmp_adata = split_adata(brm_adata) # Split into low and high BMP conditions

#### Figure 1i,j ####

# QC Analysis
low_bmp_adata_rp = preprocess(low_bmp_adata).run()

# PAGA Analysis
low_bmp_adata_rp_louvain_paga = PagaTrajectory(low_bmp_adata_rp, 'louvain').run()

# Color Coordination
color_df = color_tuple(low_bmp_adata_rp_louvain_paga, 'day_geno')
color_df['day_geno'] = ['_'.join(item.split('_')[:2]) for item in color_df.index.tolist()]
color_df = color_df.set_index('day_geno')

color_df.loc[['D4_BrmKO', 'D6_BrmKO', 'D10_BrmKO'],'new_hex'] = ['#fee6ce', '#fdae6b' ,'#e6550d']
color_df.loc[['D4_WT', 'D6_WT', 'D10_WT'],'new_hex'] = ['#deebf7', '#9ecae1', '#3182bd']

low_bmp_adata_rp_louvain_paga.uns['day_geno_colors'] = color_df['new_hex'].tolist()

#### END ####


#### Figure 3c,d ####
# QC Analysis
brm_adata_rp = preprocess(brm_adata).run()

# PAGA Analysis
brm_adata_rp_louvain_paga = PagaTrajectory(brm_adata_rp, 'leiden').run()

# Color Coordination
color_df = color_tuple(brm_adata_rp_louvain_paga, 'niceorder')

#  KO high red and low green
color_df.loc[['D4_BrmKO_highBMP4', 'D6_BrmKO_highBMP4', 'D10_BrmKO_highBMP4'],'new_hex'] = ['#fc9272', '#fb6a4a', '#de2d26']
color_df.loc[['D4_BrmKO_lowBMP4', 'D6_BrmKO_lowBMP4', 'D10_BrmKO_lowBMP4'],'new_hex'] = ['#78c679', '#31a354', '#006837']

#  WT high purple and low blue
color_df.loc[['D4_WT_highBMP4', 'D6_WT_highBMP4', 'D10_WT_highBMP4'],'new_hex'] = ['#7fcdbb', '#41b6c4', '#2c7fb8'] 
color_df.loc[['D4_WT_lowBMP4', 'D6_WT_lowBMP4', 'D10_WT_lowBMP4'],'new_hex'] = ['#df65b0', '#dd1c77','#980043']

brm_adata_rp_louvain_paga.uns['niceorder_colors'] = color_df['new_hex'].tolist()

#### END ####

# session_info
# anndata     0.7.6
# scanpy      1.8.1
# sinfo       0.3.4
# -----
# PIL                 7.1.2
# backcall            0.1.0
# cffi                1.13.2
# cloudpickle         1.3.0
# cycler              0.10.0
# cython_runtime      NA
# cytoolz             0.10.1
# dask                2.11.0
# dateutil            2.8.1
# decorator           4.4.1
# google              NA
# h5py                2.10.0
# igraph              0.8.2
# ipykernel           5.1.3
# ipython_genutils    0.2.0
# ipywidgets          7.5.1
# jedi                0.15.1
# joblib              0.15.0
# kiwisolver          1.2.0
# llvmlite            0.31.0
# louvain             0.7.0
# matplotlib          3.1.3
# more_itertools      NA
# mpl_toolkits        NA
# natsort             7.0.1
# numba               0.48.0
# numexpr             2.7.1
# numpy               1.18.4
# packaging           20.1
# pandas              1.2.4
# parso               0.5.1
# pexpect             4.7.0
# pickleshare         0.7.5
# pkg_resources       NA
# prompt_toolkit      3.0.2
# psutil              5.7.0
# ptyprocess          0.6.0
# pygments            2.5.2
# pyparsing           2.4.7
# pytz                2020.1
# scipy               1.4.1
# seaborn             0.10.1
# setuptools_scm      NA
# simplejson          3.17.2
# six                 1.13.0
# sklearn             0.22.2.post1
# statsmodels         0.11.1
# storemagic          NA
# tables              3.6.1
# tblib               1.6.0
# texttable           1.6.2
# toolz               0.10.0
# tornado             6.0.3
# traitlets           4.3.3
# wcwidth             NA
# yaml                5.3.1
# zipp                NA
# zmq                 18.1.1
# -----
# IPython             7.10.1
# jupyter_client      5.3.3
# jupyter_core        4.6.1
# jupyterlab          1.2.0
# notebook            6.0.0
# -----
# Python 3.7.6 | packaged by conda-forge | (default, Mar 23 2020, 23:03:20) [GCC 7.3.0]
# Linux-3.10.0-1160.36.2.el7.x86_64-x86_64-with-debian-buster-sid
# 72 logical CPU cores, x86_64
# -----
# Session information updated at 2021-09-13 15:43