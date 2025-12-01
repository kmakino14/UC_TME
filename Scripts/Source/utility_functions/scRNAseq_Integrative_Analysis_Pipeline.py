import sys
import os
sys.path.append(os.path.abspath('./OV_utility_functions'))
from basic_imports import *

import scanpy as sc
import infercnvpy as cnv
from scipy.sparse import vstack
import scib.metrics as sm
import anndata
import tables

from sklearn.metrics.cluster import silhouette_samples, silhouette_score
from sklearn.decomposition import PCA
pca = PCA(n_components=2)
from pingouin import linear_regression

from itertools import product
from more_itertools import unique_everseen

import anndata2ri
import rpy2.rinterface_lib.callbacks
import rpy2.robjects as ro
from scib.metrics.utils import diffusion_nn, diffusion_conn

import gseapy
import milopy
import milopy.utils
import milopy.core as milo
import milopy.plot as milopl
import scvelo as scv
import cellrank as cr

import triku as tk

from natsort import natsorted

mito_g = [i.split('.')[0] for i in list(gencode43[gencode43[0]=='chrM']['gene'].unique())]
ribo_g = gencode43.loc[[i.startswith('RPS') or i.startswith('RPL') for i in gencode43['name']], 'gene'].unique()
cyc_g = pd.read_csv('./Miscellaneous/gene_list/human_cell_cycle_signature_genes_seurat.csv')

# Integration functions

def integrate(samps, marker_genes=[], skip_cb=False, regress_out='basic', triku=False, triku_harmony=False):

    samps_to_integrate = []
    for samp in samps:
        samp_unfiltered = sc.read_h5ad('processed_anndata/%s_unfiltered.h5ad'%samp)
        samp_filtered = sc.read_h5ad('processed_anndata/%s_filtered.h5ad'%samp)
        samp_filtered_count = samp_unfiltered[samp_filtered.obs.index, ]
        samp_filtered_count.var.index.name = 'id'
        samp_filtered_count.obs = samp_filtered.obs
        samps_to_integrate.append(samp_filtered_count)

    merge = samps_to_integrate[0].concatenate(*samps_to_integrate[1:], join='outer', batch_key='sample', fill_value=0)
    merge.obs['sample'] = merge.obs['sample'].astype(int)
    for n in range(len(samps)):
        merge.obs['sample'] = merge.obs['sample'].replace(n, samps[n])
    merge.obs.index = [i.split('-')[0]+'-'+merge.obs.loc[i, 'sample'] for i in merge.obs.index]

    for i in merge.var.columns:
        if '-' in i:
            del merge.var[i]

    merge.var['mt'] = [i in mito_g for i in merge.var_names]
    merge.var['ribo'] = [i in ribo_g for i in merge.var_names]
    sc.pp.calculate_qc_metrics(merge, qc_vars=['mt', 'ribo'], percent_top=None, log1p=True, inplace=True)

    merge.write('processed_anndata/merged_unintegrated.h5ad')
    
    sc.pp.normalize_total(merge, target_sum=1e4)
    sc.pp.log1p(merge)

    if triku:

        sc.pp.filter_cells(merge, min_genes=50)
        sc.pp.filter_genes(merge, min_cells=10)

        sc.tl.pca(merge, svd_solver='arpack')

        if triku_harmony:

            merge.obsm['X_pca_uncorrected'] = merge.obsm['X_pca']
            del merge.obsm['X_pca']
            sc.external.pp.harmony_integrate(merge, key='sample', basis='X_pca_uncorrected', adjusted_basis='X_pca', max_iter_harmony=20)

        sc.pp.neighbors(merge, metric='cosine', n_neighbors=int(0.5 * len(merge) ** 0.5), n_pcs=50)

        tk.tl.triku(merge, use_raw=False)

    else:

        sc.pp.highly_variable_genes(merge, min_mean=0.0125, max_mean=3, min_disp=0.5)

    print('No. of HVG = %s'%merge.var['highly_variable'].sum())

    # Switch the indices of merge.raw from gene IDs to gene names for the convinience of downstream analysis where in most cases the gene names will be used
    merge.raw = merge
    merge.raw.var['gene_name'] = gencode43.drop_duplicates('gene').set_index('gene').reindex(index=merge.raw.var.index)['name'].tolist()
    merge_raw_indices = merge.raw.var.sort_values('total_counts', ascending=False)['gene_name'].dropna().drop_duplicates().index
    merge_raw = sc.AnnData(merge.raw[:, merge_raw_indices].X)
    merge_raw.var = merge.raw.var.loc[merge_raw_indices, ]
    merge_raw.var['gene_id'] = merge_raw.var.index
    merge_raw.var.index = merge_raw.var['gene_name']
    del merge_raw.var['gene_name']
    merge.raw = merge_raw

    sc.tl.score_genes_cell_cycle(merge, 
                                 s_genes=cyc_g.query('phase == "S"')['geneID'], 
                                 g2m_genes=cyc_g.query('phase == "G2/M"')['geneID'], 
                                 random_state=0, copy=False, use_raw=False)

    merge = merge[:, merge.var.highly_variable]

    if regress_out == 'basic':
        sc.pp.regress_out(merge, ['total_counts', 'pct_counts_mt'])
    elif regress_out == 'cell cycle':
        merge = merge[~merge.obs['S_score'].isna()]
        sc.pp.regress_out(merge, ['total_counts', 'pct_counts_mt', 'G2M_score', 'S_score'])

    sc.pp.scale(merge, max_value=10)

    sc.tl.pca(merge, svd_solver='arpack')
    merge.obsm['X_pca_uncorrected'] = merge.obsm['X_pca']
    del merge.obsm['X_pca']
    sc.pp.neighbors(merge, n_neighbors=10, n_pcs=50, use_rep='X_pca_uncorrected', key_added='uncorrected')
    sc.tl.umap(merge, neighbors_key='uncorrected')
    merge.obsm['X_umap_uncorrected'] = merge.obsm['X_umap']
    del merge.obsm['X_umap']
    sc.tl.leiden(merge, resolution=0.3, key_added='leiden_uncorrected', neighbors_key='uncorrected')

    if skip_cb:

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge, basis='X_umap_uncorrected', color=['sample', 'leiden_uncorrected', 'ploidy', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr'], \
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Sample', 'Leiden cluster (unintegrated)', 'Ploidy', 'Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR'], ncols=3)

    else:
        
        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge, basis='X_umap_uncorrected', color=['sample', 'leiden_uncorrected', 'ploidy', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr', 'cb_prob'], \
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Sample', 'Leiden cluster (unintegrated)', 'Ploidy', 'Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR', 'CellBender cell probability'], ncols=3)

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.embedding(merge, basis='X_umap_uncorrected', color=marker_genes, add_outline=False, legend_loc='right margin', \
                   legend_fontsize=8, legend_fontoutline=2, frameon=False, \
                   ncols=3, cmap=plt.cm.viridis)

    sc.external.pp.harmony_integrate(merge, key='sample', basis='X_pca_uncorrected', adjusted_basis='X_pca', max_iter_harmony=20)
    sc.pp.neighbors(merge, n_neighbors=10, n_pcs=50)
    sc.tl.umap(merge)
    sc.tl.leiden(merge, resolution=0.3)

    if skip_cb:

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge, basis='X_umap', color=['sample', 'leiden', 'ploidy', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr'], \
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Sample', 'Leiden cluster (integrated)', 'Ploidy', 'Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR'], ncols=3)
    
    else:

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge, basis='X_umap', color=['sample', 'leiden', 'ploidy', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr', 'cb_prob'], \
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Sample', 'Leiden cluster (integrated)', 'Ploidy', 'Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR', 'CellBender cell probability'], ncols=3)

    with rc_context({'figure.figsize': (5, 5)}):
        sc.pl.embedding(merge, basis='X_umap', color=marker_genes, add_outline=False, legend_loc='right margin', \
                   legend_fontsize=8, legend_fontoutline=2, frameon=False, \
                   ncols=3, cmap=plt.cm.viridis)
    
    merge.write('processed_anndata/merged_integrated.h5ad')
    return merge

def subintegrate(merge, cells_to_subset, variable_to_integrate=None, custom_genes_to_keep=[], genes_to_subset=[], regress_out='basic', triku=False, triku_harmony=False, max_iter_harmony=20):
    
    merge_sub = sc.AnnData(X=merge.raw.X, obs=merge.obs, var=merge.raw.var)
    merge_sub = merge_sub[cells_to_subset, ]
    print('Subset Data dimension: %s cells x %s genes'%merge_sub.shape)
    print('Finished.\n')
    
    print('Normalizing data...')
    merge_sub.raw = merge_sub
    if len(genes_to_subset) > 0:
        merge_sub = merge_sub[:, genes_to_subset]

    if triku:

        sc.pp.filter_cells(merge_sub, min_genes=50)
        sc.pp.filter_genes(merge_sub, min_cells=10)

        sc.tl.pca(merge_sub, svd_solver='arpack')

        if triku_harmony:

            merge_sub.obsm['X_pca_uncorrected'] = merge_sub.obsm['X_pca']
            del merge_sub.obsm['X_pca']
            sc.external.pp.harmony_integrate(merge_sub, key=variable_to_integrate, basis='X_pca_uncorrected', adjusted_basis='X_pca', max_iter_harmony=max_iter_harmony)

        sc.pp.neighbors(merge_sub, metric='cosine', n_neighbors=int(0.5 * len(merge_sub) ** 0.5), n_pcs=50)

        tk.tl.triku(merge_sub, use_raw=False)

    else:

        sc.pp.highly_variable_genes(merge_sub, min_mean=0.0125, max_mean=3, min_disp=0.5)

    print('No. of HVG = %s'%merge_sub.var['highly_variable'].sum())

    if len(custom_genes_to_keep) > 0:
        merge_sub.var.loc[custom_genes_to_keep, 'highly_variable'] = True
    merge_sub = merge_sub[:, merge_sub.var.highly_variable]

    if regress_out == 'basic':
        sc.pp.regress_out(merge_sub, ['total_counts', 'pct_counts_mt'])
    elif regress_out == 'cell cycle':
        merge_sub = merge_sub[~merge_sub.obs['S_score'].isna()]
        sc.pp.regress_out(merge_sub, ['total_counts', 'pct_counts_mt', 'G2M_score', 'S_score'])

    sc.pp.scale(merge_sub, max_value=10)
    print('Finished.\n')
    
    # Dimension reduction
    print('Dimension reduction and clustering...')
    sc.tl.pca(merge_sub, svd_solver='arpack')

    if variable_to_integrate:

        merge_sub.obsm['X_pca_uncorrected'] = merge_sub.obsm['X_pca']
        del merge_sub.obsm['X_pca']
        sc.pp.neighbors(merge_sub, n_neighbors=10, n_pcs=50, use_rep='X_pca_uncorrected', key_added='uncorrected')
        sc.tl.umap(merge_sub, neighbors_key='uncorrected')
        merge_sub.obsm['X_umap_uncorrected'] = merge_sub.obsm['X_umap']
        del merge_sub.obsm['X_umap']
        sc.tl.leiden(merge_sub, resolution=0.3, key_added='leiden_uncorrected', neighbors_key='uncorrected')
        sc.external.pp.harmony_integrate(merge_sub, key=variable_to_integrate, basis='X_pca_uncorrected', adjusted_basis='X_pca', max_iter_harmony=max_iter_harmony)

    sc.pp.neighbors(merge_sub, n_neighbors=10, n_pcs=50)
    sc.tl.umap(merge_sub)
    sc.tl.leiden(merge_sub, resolution=0.3)

    try:

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge_sub, basis='X_umap', color=['log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr'], 
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR'], ncols=2)

    except:

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(merge_sub, basis='X_umap', color=['log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo'], 
                            add_outline=False, legend_loc='right margin',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False,
                            title=['Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate'], ncols=2)


    for i in merge_sub.var.columns:
        if '-' in i:
            del merge_sub.var[i]

    print('Finished.\n')
    
    return merge_sub

# Integration quality evaluation functions

def silhouette_batch(
        adata,
        batch_key,
        group_key,
        embed,
        metric='euclidean',
        verbose=True
):
    """
    Absolute silhouette score of batch labels subsetted for each group.
    :param batch_key: batches to be compared against
    :param group_key: group labels to be subsetted by e.g. cell type
    :param embed: name of column in adata.obsm
    :param metric: see sklearn silhouette score
        default False: return average width silhouette (ASW)
    :param verbose:
    :return:
        average width silhouette ASW
        mean silhouette per group in pd.DataFrame
        Absolute silhouette scores per group label
    """
    if embed not in adata.obsm.keys():
        print(adata.obsm.keys())
        raise KeyError(f'{embed} not in obsm')

    sil_all = pd.DataFrame(columns=['group', 'silhouette_score'])

    for group in adata.obs[group_key].unique():
        adata_group = adata[adata.obs[group_key] == group]
        n_batches = adata_group.obs[batch_key].nunique()

        if (n_batches == 1) or (n_batches == adata_group.shape[0]):
            continue

        sil_per_group = silhouette_samples(
            adata_group.obsm[embed],
            adata_group.obs[batch_key],
            metric=metric
        )

        sil_all = pd.concat(
            [sil_all, 
            pd.DataFrame({
                'group': [group] * len(sil_per_group),
                'silhouette_score': sil_per_group
            })], ignore_index=True)

    sil_all = sil_all.reset_index(drop=True)

    if verbose:
        sil_means = sil_all.groupby('group').mean()
        print(f'mean silhouette per cell: {sil_means}')

    return sil_all

def kBET(adata, basis, batch, group, k0=10):
    
    anndata2ri.activate()
    ro.r("library(kBET)")
    
    groups = adata.obs[group].unique()
    df_kbet_scores = pd.DataFrame({}, columns=groups)
    
    for g in groups:
        
        adata_sub = adata[adata.obs[group]==g, ]
        #adata_tmp = diffusion_conn(adata_sub, min_k=50, copy=True)
        adata_tmp_knn = diffusion_nn(adata_sub, 10).astype(float)
        
        ro.globalenv['data_mtrx'] = adata_sub.obsm[basis]
        ro.globalenv['batch'] = adata_sub.obs[batch]

        ro.globalenv['knn_graph'] = adata_tmp_knn
        ro.globalenv['k0'] = k0
        
        ro.r(
            "batch.estimate <- kBET("
            "  data_mtrx,"
            "  batch,"
            "  knn=knn_graph,"
            "  k0=k0,"
            "  plot=FALSE,"
            "  do.pca=FALSE,"
            "  heuristic=FALSE,"
            "  adapt=FALSE,"
            ")"
        )

        try:
            kbet_scores = ro.r("batch.estimate$stats$kBET.observed")
        except rpy2.rinterface_lib.embedded.RRuntimeError:
            kbet_scores = np.nan
            
        df_kbet_scores[g] = kbet_scores

    anndata2ri.deactivate()
    return df_kbet_scores

# InferCNV functions

def inferCNV(adata, sample_col, cell_type_col, tumor_cell_type, cell_types_to_remove=[]):

    adata.raw.var['chromosome'] = gencode43.drop_duplicates('gene').set_index('gene').loc[adata.raw.var['gene_id'], 0].tolist()
    adata.raw.var['start'] = gencode43.drop_duplicates('gene').set_index('gene').loc[adata.raw.var['gene_id'], 3].tolist()
    adata.raw.var['end'] = gencode43.drop_duplicates('gene').set_index('gene').loc[adata.raw.var['gene_id'], 4].tolist()

    samples = list(adata.obs[sample_col].unique())

    for samp in samples:
        
        adata_tmp = adata[adata.obs[sample_col]==samp, ]
        adata_tmp = sc.AnnData(X = np.log((np.e**(adata_tmp.raw.X.toarray()) - 1)*100 + 1), 
                               obs = adata_tmp.obs, 
                               var = adata_tmp.raw.var)
        
        cnv.tl.infercnv(
            adata_tmp,
            reference_key=cell_type_col,
            reference_cat=list(set(adata_tmp.obs[cell_type_col]) - set(cell_types_to_remove) - set([tumor_cell_type])),
            window_size=250,
            step=1,
            exclude_chromosomes=['chrM', 'chrX', 'chrY'],
            n_jobs=32,
            inplace=True,
        )
        
        cnv_score_tmp = list(np.array(np.mean(np.abs(adata_tmp.obsm['X_cnv']), axis=1).reshape(-1))[0])
        X_cnv_tmp = adata_tmp.obsm['X_cnv']
        
        if samples.index(samp) == 0:

            cnv_score = cnv_score_tmp
            X_cnv = X_cnv_tmp

        else:

            cnv_score = cnv_score + cnv_score_tmp
            X_cnv = vstack([X_cnv, X_cnv_tmp])
            
        print(samp)
        
    adata.obs['cnv_score'] = cnv_score
    adata.obsm['X_cnv'] = X_cnv
    adata.uns['cnv'] = adata_tmp.uns['cnv']

# Mapping and scoring functions

def reference_mapping(adata, references_rds, reference_annot_labels, sample_col='sample', plot=False, use_raw=False):
    
    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.r("library(symphony)")

    if use_raw:
        ro.globalenv['query_data_mtrx'] = np.array(adata.raw.X.T)
        ro.globalenv['query_data_mtrx_colnames'] = adata.obs.index
        ro.globalenv['query_data_mtrx_rownames'] = adata.raw.var.index
    else:
        ro.globalenv['query_data_mtrx'] = adata.X.T
        ro.globalenv['query_data_mtrx_colnames'] = adata.obs.index
        ro.globalenv['query_data_mtrx_rownames'] = adata.var.index

    ro.globalenv['query_metadata'] = adata.obs[[sample_col]]
    ro.globalenv['sample_col'] = sample_col
    
    ro.r("""
    colnames(query_data_mtrx) = query_data_mtrx_colnames
    rownames(query_data_mtrx) = query_data_mtrx_rownames
    """)
    
    ref_labels = [ntpath.basename(i).split('.')[0] + '_cell_type_symphony_preds' for i in references_rds]

    for n in range(len(references_rds)):
        
        reference_rds = references_rds[n]
        reference_annot_label = reference_annot_labels[n]
        ref_label = ref_labels[n]
    
        ro.globalenv['reference_rds'] = reference_rds
        ro.globalenv['reference_annot_label'] = reference_annot_label

        ro.r("""
        reference = readRDS(reference_rds)
        """)

        ro.r("""
        query = mapQuery(query_data_mtrx, 
                     query_metadata, 
                     reference,
                     vars=c(sample_col),
                     do_normalize = FALSE,
                     do_umap = FALSE)
        """)

        ro.r("""
        query = knnPredict(query, reference, reference$meta_data[reference_annot_label][, 1], k=5)
        """)

        ro.r("""
        query_cell_types = as.vector(query$meta_data$cell_type_pred_knn)
        """)

        adata.obs[ref_label] = ro.r['query_cell_types']

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()
        
    if plot:
        
        for ref_label in ref_labels:

            adata_tmp = adata.copy()
            preds_cell_types = natsorted(adata_tmp.obs[ref_label].unique())

            for i in preds_cell_types:
                adata_tmp.obs[i] = ['Y' if j == i else 'N' for j in adata_tmp.obs[ref_label]]
                adata_tmp.obs[i] = adata_tmp.obs[i].astype(pd.CategoricalDtype(['Y', 'N'], ordered=True))

            print(ref_label.split('_cell_type_symphony_preds')[0])

            with rc_context({'figure.figsize': (4, 4), 'axes.titlesize':12}):
                sc.pl.embedding(adata_tmp, basis='X_umap', color=preds_cell_types, legend_loc=None, add_outline=False, 
                           frameon=False, title=None, palette=['darksalmon', 'lightgrey'], ncols=4)

def score_gene_sets_ucell(adata, gene_sets, plot=False, use_raw=False):

    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.r("library(UCell)")

    if use_raw:
        try:
            ro.globalenv['data_mtrx'] = adata.raw.X.T
            ro.globalenv['data_mtrx_colnames'] = adata.obs.index
            ro.globalenv['data_mtrx_rownames'] = adata.raw.var.index
        except:
            ro.globalenv['data_mtrx'] = np.array(adata.raw.X.T)
            ro.globalenv['data_mtrx_colnames'] = adata.obs.index
            ro.globalenv['data_mtrx_rownames'] = adata.raw.var.index
    else:
        ro.globalenv['data_mtrx'] = adata.X.T
        ro.globalenv['data_mtrx_colnames'] = adata.obs.index
        ro.globalenv['data_mtrx_rownames'] = adata.var.index

    ro.globalenv['gene_sets'] = [gene_sets[i] for i in list(gene_sets.keys())]

    ro.r("""
    colnames(data_mtrx) = data_mtrx_colnames
    rownames(data_mtrx) = data_mtrx_rownames
    """)

    ro.r("""
    scores = ScoreSignatures_UCell(data_mtrx, gene_sets, ncores=32)
    """)

    for i in range(ro.r['scores'].shape[1]):
        adata.obs[list(gene_sets.keys())[i] + '_gene_signature_ucell_score'] = [j[i] for j in ro.r['scores']]

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()

    if plot:

        with rc_context({'figure.figsize': (4, 4), 'axes.titlesize':12}):
            sc.pl.embedding(adata, basis='X_umap', color=sorted([i + '_gene_signature_ucell_score' for i in gene_sets.keys()]), 
                            title=sorted(gene_sets.keys()), add_outline=False, frameon=False, ncols=3, cmap=plt.cm.Purples)

# Sample-wise feature association analysis functions

def sample_pca_by_cell_frac(adata, sample_col, cell_type_col, variables, dim_reduc='pca'):

    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)
    df_frac = st.zscore(df_frac, axis=0)

    # Run PCA
    if dim_reduc == "pca":

        pca_frac = pca.fit_transform(df_frac.T)

        # Add variable columns
        df_pca_frac_annot = pd.DataFrame(pca_frac, index=df_frac.columns, columns=['PC1', 'PC2'])
        for variable in variables:
            df_pca_frac_annot[variable] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_pca_frac_annot.index, variable].tolist()

    elif dim_reduc == "umap":

        import umap
        reducer = umap.UMAP()

        pca_frac = reducer.fit_transform(df_frac.T)

        # Add variable columns
        df_pca_frac_annot = pd.DataFrame(pca_frac, index=df_frac.columns, columns=['PC1', 'PC2'])
        for variable in variables:
            df_pca_frac_annot[variable] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_pca_frac_annot.index, variable].tolist()

    variables_comb = [i for i in list(product(variables, variables)) if i[0] != i[1]]
    variables_comb = list(unique_everseen(variables_comb, key=frozenset))

    # Plot results
    if len(variables_comb) == 0:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        sns.scatterplot(data=df_pca_frac_annot, x='PC1', y='PC2', hue=variables[0], ax=ax, linewidth=0)
        axes = [ax]

    elif len(variables_comb) == 1:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6))
        sns.scatterplot(data=df_pca_frac_annot, x='PC1', y='PC2', hue=variables[0], style=variables[1], ax=ax, linewidth=0)
        axes = [ax]

    elif len(variables_comb) <= 3:
        fig, axes = plt.subplots(1, len(variables_comb), figsize=(7*len(variables_comb), 6))
        axes = np.array(axes).reshape(-1)
        for n in range(len(variables_comb)):
            ax = axes[n]
            variable1 = variables_comb[n][0]
            variable2 = variables_comb[n][1]
            sns.scatterplot(data=df_pca_frac_annot, x='PC1', y='PC2', hue=variable1, style=variable2, ax=ax, linewidth=0)

    else:
        nrows = len(variables_comb)//3+1 if len(variables_comb)%3>0 else len(variables_comb)//3
        fig, axes = plt.subplots(nrows, 3, figsize=(21, 6*nrows))
        axes = np.array(axes).reshape(-1)
        for n in range(len(variables_comb)):
            ax = axes[n]
            variable1 = variables_comb[n][0]
            variable2 = variables_comb[n][1]
            sns.scatterplot(data=df_pca_frac_annot, x='PC1', y='PC2', hue=variable1, style=variable2, ax=ax, linewidth=0)

    for ax in axes:

        sns.despine(right=True, top=True, ax=ax)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Add legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        if len(variables_comb) == 0:
            legend = ax.legend(handles, labels, title=variables[0], loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}, frameon=False)
            legend.get_title().set_size(12)
        else:
            legend = ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 8}, frameon=False)

def cell_frac_feature_imp(adata, sample_col, cell_type_col, variables, plot=True, linewidth=3, figsize=(5, 5), colormap=plt.cm.magma):

    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)
    df_frac = st.zscore(df_frac, axis=0)

    X = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac.columns, variables]
    for i in X.columns:
        X[i] = X[i].tolist()
    X = pd.get_dummies(X, drop_first=True)

    df_relimp = pd.DataFrame({}, index=df_frac.index, columns=X.columns)

    for cell_type in df_frac.index:

        y = df_frac.loc[cell_type, ]
        res = linear_regression(X, y, remove_na=True, relimp=True)
        df_relimp = df_relimp.loc[:, res.names[1:]]
        df_relimp.loc[cell_type, ] = res.iloc[1:, ]['relimp'].tolist()

    df_relimp = df_relimp[sorted(df_relimp.columns)]

    if plot:

        cm = sns.clustermap(df_relimp.astype(float), col_cluster=False, annot=False, figsize=figsize,
                            linewidth=linewidth, vmax=min(0.5, df_relimp.max().max()), cmap=colormap)
        cm.ax_cbar.set_ylabel('Feature importance', fontsize=14)
        cm.ax_cbar.yaxis.set_label_position("left")
        cm.ax_cbar.yaxis.tick_left()
        cm.ax_heatmap.set_xticklabels(df_relimp.columns, rotation=45, ha='right')

    else:

        return df_relimp

# Differential abundance testing functions

def filter_adata_by_cell_count(adata, sample_col, threshold):

    filtered_samps = adata.obs.groupby(sample_col).count().query('total_counts >= @threshold').index.tolist()
    adata_filtered = adata[adata.obs.query('sample in @filtered_samps').index]

    return adata_filtered

def run_milo(adata, sample_col, cell_type_col, variable_to_test, covariates, level_to_be_contrasted=None, samples_to_subset=None, 
             basis='X_umap', prop=0.3, alpha=0.25, min_logFC=1, cmap=plt.cm.RdBu_r, plot_edges=False, save=None):

    #if 'nhood_adata' not in list(adata.uns.keys()):
    #    milo.make_nhoods(adata, prop=prop)
    #    milo.count_nhoods(adata, sample_col="sample")

    #if samples_to_subset:
    #    adata_ = adata[adata.obs[sample_col].isin(samples_to_subset)].copy()
    #else:
    #    adata_ = adata.copy()

    #sc.tl.pca(adata_, svd_solver='arpack')
    #sc.pp.neighbors(adata_, n_neighbors=10, n_pcs=50)

    adata_ = adata.copy()

    variable_to_test_orig = variable_to_test
    if ' ' in variable_to_test:
        variable_to_test = variable_to_test.replace(' ', '_')
        adata_.obs[variable_to_test] = [i.replace(' ', '_') for i in adata_.obs[variable_to_test_orig]]
        adata_.obs[variable_to_test] = adata_.obs[variable_to_test].astype('category')

    level_to_be_contrasted_orig = level_to_be_contrasted
    if ' ' in level_to_be_contrasted:
        level_to_be_contrasted_orig = level_to_be_contrasted
        level_to_be_contrasted = level_to_be_contrasted.replace(' ', '_')

    if level_to_be_contrasted:
        adata_.obs[variable_to_test] = adata_.obs[variable_to_test].cat.reorder_categories(list(set(adata_.obs[variable_to_test].dropna().unique()) - set([level_to_be_contrasted]))+[level_to_be_contrasted])
        title = '%s log2FC: %s vs. the rest'%(variable_to_test_orig, level_to_be_contrasted_orig)
    else:
        title = '%s correlation'%variable_to_test

    milo.make_nhoods(adata_, prop=prop)
    milo.count_nhoods(adata_, sample_col=sample_col)

    milo.DA_nhoods(adata=adata_, 
                   design="~%s+%s"%(variable_to_test, '+'.join(covariates)) if covariates else "~%s"%variable_to_test, 
                   subset_samples=set(adata_.obs[sample_col].unique()) & set(samples_to_subset) if samples_to_subset else set(adata_.obs[sample_col].unique()))

    milopy.utils.build_nhood_graph(adata_, basis=basis)
    milopy.utils.annotate_nhoods(adata_, anno_col=cell_type_col)
    
    def plot_nhood_graph(
        adata,
        alpha = 0.1,
        min_logFC = 0,
        min_size = 10,
        plot_edges = False,
        title = "DA log-Fold Change",
        **kwargs
    ):
        nhood_adata = adata.uns["nhood_adata"].copy()

        if "Nhood_size" not in nhood_adata.obs.columns:
            raise KeyError(
                'Cannot find "Nhood_size" column in adata.uns["nhood_adata"].obs -- \
                    please run milopy.utils.build_nhood_graph(adata)'
            )

        nhood_adata.obs["graph_color"] = nhood_adata.obs["logFC"]
        nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"]
                            > alpha, "graph_color"] = np.nan
        nhood_adata.obs["abs_logFC"] = abs(nhood_adata.obs["logFC"])
        nhood_adata.obs.loc[nhood_adata.obs["abs_logFC"]
                            < min_logFC, "graph_color"] = np.nan

        # Plotting order - extreme logFC on top
        nhood_adata.obs.loc[nhood_adata.obs["graph_color"].isna(),
                            "abs_logFC"] = np.nan
        ordered = nhood_adata.obs.sort_values(
            'abs_logFC', na_position='first').index
        nhood_adata = nhood_adata[ordered]

        vmax = np.max([nhood_adata.obs["graph_color"].max(),
                      abs(nhood_adata.obs["graph_color"].min())])
        vmin = - vmax

        sc.pl.embedding(nhood_adata, "X_milo_graph",
                        color="graph_color", cmap=cmap,
                        size=adata.uns["nhood_adata"].obs["Nhood_size"]*min_size,
                        edges=plot_edges, neighbors_key="nhood",
                        # edge_width =
                        sort_order=False,
                        frameon=False,
                        vmax=vmax, vmin=vmin,
                        title=title,
                        **kwargs
                        )

    with rc_context({'figure.figsize': (5, 5)}):
        plot_nhood_graph(adata_, 
                        alpha=alpha,
                        min_logFC=min_logFC, 
                        min_size=1,
                        title=title,
                        plot_edges=plot_edges,
                        show=False)
    if save:
        plt.savefig(save+'_milo_nhood_graph.pdf', transparent=True)

    plt.figure(figsize=(6, 8))
    milopl.plot_DA_beeswarm(
        adata_, 
        anno_col = 'nhood_annotation', 
        alpha=alpha,
        subset_nhoods=None
    )
    sns.despine(right=True, top=True)
    plt.xlabel(title)
    plt.ylabel('')
    if save:
        plt.savefig(save+'_milo_DA_beeswarm.pdf', transparent=True)
    plt.show()

# Differential gene expression testing functions

def singlecell_de(adata, sample_col, 
                  variable_to_test, level_to_be_contrasted, level_to_contrast, 
                  covariate=None, covariate_subset=None,
                  log2fc=1, pval=0.05):

    if covariate:

        adata_ = adata[adata.obs[covariate]==covariate_subset]

    else:

        adata_ = adata.copy()

    sc.tl.rank_genes_groups(adata_, variable_to_test, groups=[level_to_be_contrasted], reference=level_to_contrast, method="t-test")

    df_de_res = sc.get.rank_genes_groups_df(adata_, group=level_to_be_contrasted)

    df_de_res['pvalue'] = df_de_res['pvals']
    del df_de_res['pvals']

    df_de_res['log2FoldChange'] = df_de_res['logfoldchanges']
    del df_de_res['logfoldchanges']

    df_de_res['padj'] = df_de_res['pvals_adj']
    del df_de_res['pvals_adj']

    df_de_res['-log10pvals'] = -np.log10(df_de_res['pvalue'])
    df_de_res['significant'] = ['True' if abs(df_de_res.loc[i, 'log2FoldChange'])>log2fc and df_de_res.loc[i, 'padj']<pval else 'False' for i in df_de_res.index]
    df_de_res = df_de_res.rename(columns={'names':'name'})

    return df_de_res

def pseudobulk_de(adata, count_adata, sample_col, 
                  variable_to_test, level_to_be_contrasted, level_to_contrast, 
                  covariate=None, sample_dict=None, sample_pairing_dict=None,
                  log2fc=1, pval=0.05, lfcshrink=True):
    
    if not sample_dict:

        adata_sub = adata[(adata.obs[variable_to_test].isin([level_to_be_contrasted, level_to_contrast]))]

    else:

        adata_sub = adata[adata.obs[sample_col].isin([a for b in list(sample_dict.values()) for a in b])]
        adata_sub.obs[variable_to_test] = np.nan
        for i in list(sample_dict.keys()):
            adata_sub.obs.loc[adata_sub.obs[sample_col].isin(sample_dict[i]), variable_to_test] = i
        adata_sub.obs[variable_to_test] = adata_sub.obs[variable_to_test].astype('category')

    if sample_pairing_dict:

        adata_sub.obs[covariate] = adata_sub.obs[sample_col].map(sample_pairing_dict)
        paired_samps = adata_sub.obs.drop_duplicates(sample_col).groupby(covariate).count().query('total_counts == 2').index
        adata_sub = adata_sub[adata_sub.obs[covariate].isin(paired_samps)]

    adata_sub.obs[variable_to_test] = adata_sub.obs[variable_to_test].cat.reorder_categories([level_to_contrast, level_to_be_contrasted])

    df_count = pd.DataFrame(count_adata[adata_sub.obs.index].X.toarray(), \
                                   index=adata_sub.obs.index, \
                                   columns=count_adata.var.index)
    adata_sub.obs['Identifier'] = (adata_sub.obs[sample_col].astype(str) + \
                                    '-' + \
                                    adata_sub.obs[variable_to_test].astype(str)).tolist()
    df_count['Identifier'] = adata_sub.obs['Identifier'].tolist()
    df_count_bulk = df_count.groupby('Identifier').sum().T
    df_cpm_bulk = df_count_bulk.divide(df_count_bulk.sum(axis=0), axis=1)*1e6
    
    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.r("library(DESeq2)")
    ro.r("library(BiocParallel)")
    
    ro.globalenv['data_mtrx'] = df_count_bulk.to_numpy()
    ro.globalenv['data_mtrx_colnames'] = df_count_bulk.columns
    ro.globalenv['data_mtrx_rownames'] = df_count_bulk.index

    ro.r("""
    colnames(data_mtrx) = data_mtrx_colnames
    rownames(data_mtrx) = data_mtrx_rownames
    """)
            
    if covariate:
        
        ro.globalenv['meta'] = adata_sub.obs.drop_duplicates('Identifier').set_index('Identifier')[[sample_col, variable_to_test, covariate]]
        ro.r("""
        data_mtrx = data_mtrx[, as.vector(rownames(meta))]
        dds = DESeqDataSetFromMatrix(countData = data_mtrx, colData = meta, design = ~ %s + %s)
        """%(covariate.replace('-', '.'), variable_to_test.replace('-', '.')))
        
    else:
        
        ro.globalenv['meta'] = adata_sub.obs.drop_duplicates('Identifier').set_index('Identifier')[[sample_col, variable_to_test]]
        ro.r("""
        data_mtrx = data_mtrx[, as.vector(rownames(meta))]
        dds = DESeqDataSetFromMatrix(countData = data_mtrx, colData = meta, design = ~ %s)
        """%variable_to_test.replace('-', '.'))
    
    if lfcshrink:

        ro.r("""
        dds <- DESeq(dds)
        res = lfcShrink(dds, coef="%s_%s_vs_%s", type="apeglm")
        """%(variable_to_test.replace('-', '.'), level_to_be_contrasted.replace('-', '.'), level_to_contrast.replace('-', '.')))

        deseq_res = pd.DataFrame(ro.r('res@listData'), 
                                 index=['baseMean', 'log2FoldChange', 'lfcSE', 'pvalue', 'padj'], 
                                 columns=ro.r('res@rownames')).T

    else:

        ro.r("""
        dds <- DESeq(dds)
        res = results(dds, name="%s_%s_vs_%s")
        """%(variable_to_test.replace('-', '.'), level_to_be_contrasted.replace('-', '.'), level_to_contrast.replace('-', '.')))
    
        deseq_res = pd.DataFrame(ro.r('res@listData'), 
                                 index=['baseMean', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue', 'padj'], 
                                 columns=ro.r('res@rownames')).T
    
    if deseq_res.index[0].startswith('ENSG'):
        deseq_res['name'] = gencode43.drop_duplicates('gene').set_index('gene').reindex(index=deseq_res.index)['name'].tolist()
    else:
        deseq_res['name'] = deseq_res.index.tolist()
    
    deseq_res['-log10pvals'] = -np.log10(deseq_res['pvalue'])

    padj_thresh = min(deseq_res['padj'].dropna(), key=lambda x:abs(x - pval))
    pval_thresh = deseq_res.query('padj == @padj_thresh')['pvalue'][0]
    deseq_res['significant'] = ['True' if (abs(deseq_res.loc[i, 'log2FoldChange'])>log2fc) 
                                and (deseq_res.loc[i, 'pvalue']<pval_thresh) 
                                else 'False' for i in deseq_res.index]
    
    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()
    
    return df_cpm_bulk, deseq_res

def run_pathway_enrichment(df_de_res, terms=['GO_Biological_Process_2021', 'WikiPathways_2019_Human'], log2fc=1, pval=0.05):

    padj_thresh = min(df_de_res['padj'].dropna(), key=lambda x:abs(x - pval))
    pval_thresh = df_de_res.query('padj == @padj_thresh')['pvalue'].iloc[0]

    df_de_res = df_de_res.sort_values('pvalue').drop_duplicates('name')

    if 'scores' in df_de_res.columns.tolist():
        gene_rank = df_de_res.query('log2FoldChange > @log2fc or log2FoldChange < -@log2fc').query('pvalue < @pval_thresh')[['name', 'scores']]
    else:
        gene_rank = df_de_res.query('log2FoldChange > @log2fc or log2FoldChange < -@log2fc').query('pvalue < @pval_thresh')[['name', 'log2FoldChange']]

    enr_res_dict = {}
    gsea_res_dict = {}

    for term in terms:

        enr_res = gseapy.enrichr(gene_list=conv_symbol(df_de_res.query('log2FoldChange > @log2fc and pvalue < @pval_thresh')['name'].tolist()),
                                 organism='Human',
                                 gene_sets=term,
                                 description='pathway',
                                 cutoff = 0.5)
        enr_res_dict[term] = enr_res

        gsea_res = gseapy.prerank(rnk=gene_rank, gene_sets=term, processes=32)
        gsea_res_dict[term] = gsea_res

    return enr_res_dict, gsea_res_dict

def run_slingshot(adata, basis, cell_type_col, start_clust, end_clust, plot=True):

    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.r("library(slingshot)")

    ro.globalenv['data_mtrx'] = adata.obsm[basis]
    ro.globalenv['clust_labels'] = adata.obs[cell_type_col]
    ro.globalenv['start_clust'] = start_clust
    ro.globalenv['end_clust'] = end_clust

    ro.r("""
    traj_sling = slingshot(data=data_mtrx, \
                           clusterLabels=clust_labels, \
                           start.clus=start_clust, \
                           end.clus=end_clust, \
                           use.median=TRUE)
    """)

    dict_traj_curves = {}

    lineages = ['Lineage%s'%n for n in range(1, len(ro.r('traj_sling@metadata$curves'))+1)]

    for lineage in lineages:
        start = ro.r('traj_sling@metadata$lineages$%s'%lineage)[0].replace(' ', '-')
        end = ro.r('traj_sling@metadata$lineages$%s'%lineage)[-1].replace(' ', '-')
        dict_traj_curves[lineage] = pd.DataFrame(ro.r('traj_sling@metadata$curves$%s$s'%lineage), columns=['x', 'y'])
        adata.obs['traj_sling_time_%s_%s_to_%s'%(lineage, start, end)] = [i[lineages.index(lineage)] for i in ro.r('traj_sling@assays@data@listData$pseudotime')]

    if plot:

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))

        with rc_context({'figure.figsize': (5, 5)}):
            sc.pl.embedding(adata, basis=basis, color=[cell_type_col], \
                            add_outline=False, legend_loc='on data',
                            legend_fontsize=8, legend_fontoutline=2, frameon=False, show=False, ax=ax)

        for lineage in lineages:
            
            ax.scatter(x=dict_traj_curves[lineage]['x'], \
                       y=dict_traj_curves[lineage]['y'], \
                       s=1, linewidth=1, c='k')

        plt.show()

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()

    return dict_traj_curves

def run_cellchat(adatas, labels, cell_type_col, count_adata, cc_obj_path, cc_obj_label, downsample_n_cell=5000):

    adatas_ = []

    for n in range(len(adatas)):

        adata = adatas[n].copy()
        adata.obs[cell_type_col[n]] = [(i + ' ' + labels[n]).strip() for i in adata.obs[cell_type_col[n]]]
        adatas_.append(adata)

    cell_indices = list(set([a for b in [adata.obs.index.tolist() for adata in adatas] for a in b]))

    count_adata_ = count_adata[cell_indices, ]
    count_adata_.obs['cell_type'] = np.nan
    count_adata_.obs['cell_state'] = np.nan

    for n in range(len(adatas_)):

        adata = adatas_[n]
        count_adata_.obs.loc[adata.obs.index, 'cell_type'] = cc_obj_label.split('_')[n]
        count_adata_.obs.loc[adata.obs.index, 'cell_state'] = adata.obs[cell_type_col[n]].tolist()

    count_adata_ = count_adata_[np.random.choice(range(len(count_adata_)), downsample_n_cell, replace=False), ]

    count_adata_.obs['cell_state'] = count_adata_.obs['cell_state'].astype('category')
    count_adata_.obs['cell_state'] = count_adata_.obs['cell_state'].cat.reorder_categories(list(count_adata_.obs.sort_values('cell_type')['cell_state'].unique()))

    sc.pp.normalize_total(count_adata_, target_sum=1e4)
    sc.pp.log1p(count_adata_)
    sc.pp.highly_variable_genes(count_adata_, min_mean=0.0125, max_mean=3, min_disp=0.5)
    count_adata_ = count_adata_[:, count_adata_.var.highly_variable]

    count_adata_.var.index = gencode43.drop_duplicates('gene').set_index('gene').reindex(index=count_adata_.var.index)['name'].tolist()
    count_adata_ = count_adata_[:, ~count_adata_.var.index.duplicated()]
    count_adata_.var.index = gencode43.drop_duplicates('gene').set_index('gene').reindex(index=count_adata_.var.index)['name'].tolist()

    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.r("library(SeuratDisk)")
    ro.r("library(CellChat)")

    ro.globalenv['data_mtrx'] = count_adata_.X.T
    ro.globalenv['data_mtrx_colnames'] = count_adata_.obs.index
    ro.globalenv['data_mtrx_rownames'] = count_adata_.var.index

    ro.globalenv['meta_data'] = count_adata_.obs[['cell_state']]
    ro.globalenv['cc_obj_path'] = cc_obj_path
    ro.globalenv['cc_obj_label'] = cc_obj_label

    ro.r("""
    colnames(data_mtrx) = data_mtrx_colnames
    rownames(data_mtrx) = data_mtrx_rownames
    """)

    ro.r("""

        cellchat <- createCellChat(object = data_mtrx, 
                                   meta = meta_data, group.by = "cell_state")
                                   
        CellChatDB <- CellChatDB.human
        CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
        cellchat@DB <- CellChatDB.use

        cellchat <- subsetData(cellchat)
        #future::plan("multiprocess", workers = 32)

        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)

        cellchat <- computeCommunProb(cellchat)
        cellchat <- filterCommunication(cellchat, min.cells = 10)

        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)

        saveRDS(cellchat, file = paste0(cc_obj_path, sprintf("/cc_obj_%s.rds", cc_obj_label)))

        """)

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()

# Network analysis

from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell

RESOURCES_FOLDER="/priv18data1/hliang1_group/yluo5/Miscellaneous/scenic_auxiliaries"
DATABASE_FOLDER = RESOURCES_FOLDER+"/feather/"

DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38_*_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = os.path.join(RESOURCES_FOLDER, "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl")
MM_TFS_FNAME = os.path.join(RESOURCES_FOLDER, 'allTFs_hg38.txt')

tf_names = load_tf_names(MM_TFS_FNAME)
db_fnames = glob.glob(DATABASES_GLOB)
def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]

def run_scenic(adata, adata_count, outdir, nes_threshold=3):

    db1 = pd.read_feather(DATABASE_FOLDER+'hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
    db2 = pd.read_feather(DATABASE_FOLDER+'hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather')
    motif = pd.read_table(MOTIF_ANNOTATIONS_FNAME)
    tfs = pd.read_table(MM_TFS_FNAME, index_col=0, header=None)

    scenic_genes = list(set(db1.columns)|set(db2.columns)|set(motif['gene_name'])|set(tfs.index))

    MODULES_FNAME = os.path.join(outdir, 'modules.p')
    MOTIFS_FNAME = os.path.join(outdir, "motifs.csv")
    REGULONS_FNAME = os.path.join(outdir, "regulons.p")
    ADATA_FNAME = os.path.join(outdir, "auc.h5ad")

    adata_ = adata_count[list(set(adata.obs.index)&set(adata_count.obs.index)), ].copy()
    sc.pp.filter_genes(adata_, min_cells=5)
    adata_.var.index = gencode43.drop_duplicates('gene').set_index('gene').reindex(index=adata_.var.index)['name'].tolist()
    adata_ = adata_[:, ~adata_.var.index.duplicated()]
    adata_ = adata_[:, list(set(scenic_genes)&set(adata_.var.index))]

    df_adata = pd.DataFrame(adata_.X.toarray(), 
                            index=adata_.obs.index,
                            columns=adata_.var.index)

    # Infer co-expression TF modules
    adjacencies = grnboost2(expression_data=adata_.X, gene_names=adata_.var.index, tf_names=tf_names, verbose=True)
    modules = list(modules_from_adjacencies(adjacencies, df_adata))
    with open(MODULES_FNAME, 'wb') as f:
        pickle.dump(modules, f)
        f.close()

    # Calculate a list of enriched motifs and the corresponding target genes for all modules.
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME, nes_threshold=nes_threshold, num_workers=40, client_or_address='dask_multiprocessing')
    df.to_csv(MOTIFS_FNAME)

    # Create regulons from this table of enriched motifs.
    regulons = df2regulons(df)
    with open(REGULONS_FNAME, "wb") as f:
        pickle.dump(regulons, f)

    auc_mtx = aucell(df_adata, regulons, num_workers=40)
    
    auc_adata = anndata.AnnData(X=auc_mtx.to_numpy(),
                                obs=pd.DataFrame(auc_mtx.index),
                                var=pd.DataFrame(auc_mtx.columns))

    auc_adata.var.index = auc_adata.var['Regulon']
    auc_adata.obs.index = auc_adata.obs['Cell']
    auc_adata.write(ADATA_FNAME)

    return auc_adata