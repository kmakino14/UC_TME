import sys
import os
sys.path.append(os.path.abspath('./OV_utility_functions'))
from basic_imports import *

import scanpy as sc
import anndata
import tables

from scipy.io import mmwrite
from scipy.io import mmread
from scipy.sparse import csr_matrix

import rpy2
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
pandas2ri.activate()

droputil = importr('DropletUtils')
copykat = importr('copykat')
soupx = importr('SoupX')
Matrix = importr('Matrix')

counts = ro.r['counts']
dim = ro.r['dim']
set_seed = ro.r['set.seed']
metadata = ro.r['metadata']
as_matrix = ro.r['as.matrix']
as_character = ro.r['as.character']

mito_g = [i.split('.')[0] for i in list(gencode43[gencode43[0]=='chrM']['gene'].unique())]
ribo_g = gencode43.loc[[i.startswith('RPS') or i.startswith('RPL') for i in gencode43['name']], 'gene'].unique()
cyc_g = pd.read_csv('./Miscellaneous/gene_list/human_cell_cycle_signature_genes_seurat.csv')

class HaltException(Exception): pass

def basic_qc_plot(axes, adata, ed_low):
    
    # Saturation plot
    x=np.log10(adata.obs['total_counts']+1)
    y=np.log10(adata.obs['n_genes_by_counts']+1)

    sns.scatterplot(x=x, y=y, s=0.5, alpha=0.7, linewidth=0, ax=axes[0, 0], rasterized=True)
    axes[0, 0].set_xlabel('UMI count per cell (log10)')
    axes[0, 0].set_ylabel('Gene count per cell (log10)')
    
    # Knee plot: UMI count
    x=np.log10(range(1, adata.X.shape[0]+1))
    y=sorted(np.log10(adata.obs['total_counts']+1))[::-1]
    knee = adata.uns['knee']
    inflec = adata.uns['inflec']

    sns.scatterplot(x=x, y=y, s=1, linewidth=0, ax=axes[0, 1], rasterized=True)
    axes[0, 1].axhline(y=np.log10(knee), ls='--', lw=1, c='salmon', label='Knee point (%s)'%int(knee))
    axes[0, 1].axhline(y=np.log10(inflec), ls='--', lw=1, c='green', label='Inflection point (%s)'%int(inflec))
    axes[0, 1].set_xlabel('UMI count rank (log10)')
    axes[0, 1].set_ylabel('UMI count per cell (log10)')
    axes[0, 1].legend(loc=3)
    
    # Knee plot: gene count
    x=np.log10(range(1, adata.X.shape[0]+1))
    y=sorted(np.log10(adata.obs['n_genes_by_counts']+1))[::-1]

    sns.scatterplot(x=x, y=y, s=1, linewidth=0, ax=axes[0, 2], rasterized=True)
    #plt.axhline(y=np.log10(1000), ls='--', lw=1, c='k')
    axes[0, 2].set_xlabel('Gene count rank (log10)')
    axes[0, 2].set_ylabel('Gene count per cell (log10)')
    
    # Mito reads fraction
    sns.kdeplot(adata.obs['pct_counts_mt'], ax=axes[1, 0])
    axes[1, 0].set_xlabel('Porportion of mtDNA UMI counts')
    
    # miQC plot
    x=adata.obs['n_genes_by_counts']
    y=adata.obs['pct_counts_mt']

    sns.scatterplot(x=x, y=y, s=0.5, alpha=0.5, linewidth=0, ax=axes[1, 1], rasterized=True)
    axes[1, 1].set_xlabel('Gene count per cell')
    axes[1, 1].set_ylabel('Porportion of mtDNA UMI counts')
    
    # Print empty droplets
    x=np.log10(adata.obs['total_counts']+1)
    y=-adata.obs['ed_logprob']

    sns.scatterplot(x=x, y=y, s=0.5, hue=adata.obs['ed_class'], hue_order=['Empty droplet', 'Cell'], linewidth=0, ax=axes[1, 2], rasterized=True)
    axes[1, 2].axvline(x=np.log10(knee), ls='--', lw=1, c='k', label='Knee point (%s)'%int(knee))
    axes[1, 2].set_xlabel('UMI count per cell (log10)')
    axes[1, 2].set_ylabel('EmptyDrops -logProb')
    axes[1, 2].legend(title='')

    # Print input for CellBender
    print('No. of barcodes above the knee point (nUMI > %s): %s'%(knee, adata.obs.query('total_counts > %s'%knee).shape[0]))
    print('No. of barcodes above the lower point (nUMI > %s): %s'%(ed_low, adata.obs.query('total_counts > @ed_low').shape[0]))
    
    return adata

def preprocessing(adata_filtered):
    
    # Normalization
    print('Normalizing data...')
    sc.pp.normalize_total(adata_filtered, target_sum=1e4)
    sc.pp.log1p(adata_filtered)
    sc.pp.highly_variable_genes(adata_filtered, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata_filtered.raw = adata_filtered
    adata_filtered = adata_filtered[:, adata_filtered.var.highly_variable]
    sc.pp.regress_out(adata_filtered, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adata_filtered, max_value=10)
    print('Finished.\n')
    
    # Dimension reduction
    print('Dimension reduction and clustering...')
    sc.tl.pca(adata_filtered, svd_solver='arpack')
    sc.pp.neighbors(adata_filtered, n_neighbors=10, n_pcs=50)
    sc.tl.umap(adata_filtered)
    sc.tl.leiden(adata_filtered, resolution=0.5)
    print('Finished.\n')
    
    return adata_filtered

def advanced_qc_plot(adata_filtered, marker_genes, filter='EmptyDrops'):
    
    if filter == 'CellBender':
        with rc_context({'figure.figsize': (3, 3)}):
            fig_meta = sc.pl.embedding(adata_filtered, basis='X_umap', color=['leiden', 'log1p_total_counts', 'pct_counts_mt', 'doublet_score', 'ed_fdr', 'cb_prob'], \
                                        add_outline=False, legend_loc='on data',
                                        legend_fontsize=8, legend_fontoutline=2, frameon=False,
                                        title=['Leiden cluster', 'Total UMI count', 'Mitochondrial UMI rate', 'Doublet score', 'EmptyDrops cell FDR', 'CellBender cell probability'], ncols=3, return_fig=True)
    else:
        with rc_context({'figure.figsize': (3, 3)}):
            fig_meta = sc.pl.embedding(adata_filtered, basis='X_umap', color=['leiden', 'log1p_total_counts', 'pct_counts_mt', 'pct_counts_ribo', 'doublet_score', 'ed_fdr'], \
                                        add_outline=False, legend_loc='on data',
                                        legend_fontsize=8, legend_fontoutline=2, frameon=False,
                                        title=['Leiden cluster', 'Total UMI count', 'Mitochondrial UMI rate', 'Ribosomal UMI rate', 'Doublet score', 'EmptyDrops cell FDR'], ncols=3, return_fig=True)
    
    marker_genes = [i for i in marker_genes if i in adata_filtered.raw.var['gene_name'].tolist()]
    ids_to_plot = gencode43.drop_duplicates('name').set_index('name').loc[marker_genes, 'gene']

    with rc_context({'figure.figsize': (3, 3)}):
        fig_genes = sc.pl.embedding(adata_filtered, basis='X_umap', color=ids_to_plot, add_outline=False, legend_loc='right margin', \
                                     legend_fontsize=8, legend_fontoutline=2, frameon=False, 
                                     ncols=3, title=marker_genes, cmap=plt.cm.viridis, return_fig=True)

    return fig_meta, fig_genes

def run_pipeline_one(sample, marker_genes=None, quantifier='CellRanger', plot_only=False, save=None, filter='EmptyDrops', ed_fdr=0.01, ed_low=500, inflection=False):
    
    # Reads unfiltered AnnData
    print('Reading unfiltered data of %s from %s'%(sample, quantifier))
    if quantifier == 'CellRanger':
        try:
            adata = sc.read_10x_h5('cellranger_out/%s/outs/raw_feature_bc_matrix.h5'%sample)
        except:
            adata = sc.read_10x_mtx('cellranger_out/%s/outs/raw_feature_bc_matrix'%sample)
        adata.var['gene_name'] = adata.var.index
        adata.var.index = adata.var['gene_ids']
        del adata.var['gene_ids']
        adata.obs.index = [i.split('-')[0] for i in adata.obs.index]
    elif quantifier == 'KB':
        adata = sc.read_h5ad('kb_out/%s/counts_unfiltered/adata.h5ad'%sample)
    elif quantifier == "CeleScope":
        csdir = 'celescope_out/%s/outs/raw_feature_bc_matrix/'%sample
        if '.gz' in glob.glob(csdir+'/*')[0]:
            subprocess.call('gunzip '+csdir+'/*', shell=True)
        adata = sc.read_mtx(csdir+'matrix.mtx')
        adata = adata.T
        adata.var = pd.read_table(csdir+'genes.tsv', header=None, index_col=0, names=['id', 'name'])
        adata.var['gene_name'] = adata.var['name']
        del adata.var['name']
        adata.obs = pd.read_table(csdir+'barcodes.tsv', header=None, index_col=0, names=['barcode'])
    else:
        raise HaltException('Quantifier %s is not supported!'%quantifier)
    adata.uns['sample'] = sample
    print('Finished.\n')
    
    # Compute basic QC metrics
    adata.var.index = [i.split('.')[0] for i in adata.var.index]
    adata = adata[:, ~adata.var.index.duplicated()]
    adata = adata[:, list(set(adata.var.index)&set(gencode43['gene']))]
    adata.var['mt'] = [i in mito_g for i in adata.var_names]
    adata.var['ribo'] = [i in ribo_g for i in adata.var_names]
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo'], percent_top=None, log1p=True, inplace=True)
    
    # Get knee and inflection points
    print('Running EmptyDrops...')
    if quantifier == 'CellRanger':

        crdir = 'cellranger_out/%s/outs/'%sample
        mmwrite(crdir+'/matrix.mtx', csr_matrix(adata.X.T), field='integer')
        adata.obs.reset_index().iloc[:, 0].to_csv(crdir+'/barcodes.tsv', header=None, sep='\t', index=False)
        gencode43.drop_duplicates('gene').set_index('gene').loc[adata.var.index, 'name'].to_csv(crdir+'/genes.tsv', header=None, sep='\t')
        tmp = droputil.read10xCounts(crdir)

        #if '.gz' in glob.glob('cellranger_out/%s/outs/raw_feature_bc_matrix/*'%sample)[0]:
        #    for i in glob.glob('cellranger_out/%s/outs/raw_feature_bc_matrix/*gz'%sample):
        #        subprocess.call('gunzip -c ' + i + '>%s'%i.replace('.gz', ''), shell=True)
        #tmp = droputil.read10xCounts('cellranger_out/%s/outs/raw_feature_bc_matrix'%sample)
    
    elif quantifier  == 'KB':
        kbdir = 'kb_out/%s/counts_unfiltered/'%sample
        mmwrite(kbdir+'/matrix.mtx', csr_matrix(adata.X.T), field='integer')
        adata.obs.reset_index().iloc[:, 0].to_csv(kbdir+'/barcodes.tsv', header=None, sep='\t', index=False)
        gencode43.drop_duplicates('gene').set_index('gene').loc[adata.var.index, 'name'].to_csv(kbdir+'/genes.tsv', header=None, sep='\t')
        tmp = droputil.read10xCounts(kbdir)
    elif quantifier == "CeleScope":
        tmp = droputil.read10xCounts('celescope_out/%s/outs/raw_feature_bc_matrix'%sample)
    else:
        raise HaltException('Quantifier %s is not supported!'%quantifier)
        
    set_seed(0)
    tmp_br_out = droputil.barcodeRanks(counts(tmp), lower=ed_low)
    adata.uns['knee'] = int(metadata(tmp_br_out)[0][0])
    adata.uns['inflec'] = int(metadata(tmp_br_out)[1][0])

    # Call empty droplets
    set_seed(10)
    if inflection:
        ed_low = adata.uns['inflec']
    else:
        ed_low = ed_low
    tmp_e_out = droputil.emptyDrops(counts(tmp), lower=ed_low)
    adata.obs['ed_logprob'] = pd.DataFrame(as_matrix(tmp_e_out)).iloc[:, 1].tolist()
    adata.obs['ed_fdr'] = pd.DataFrame(as_matrix(tmp_e_out)).iloc[:, -1].tolist()
    adata.obs['ed_class'] = ['Cell' if i < ed_fdr else 'Empty droplet' for i in adata.obs['ed_fdr']]
    print('Finished.\n')
    
    if plot_only:

        fig, axes = plt.subplots(2, 3, figsize=(16, 8))
        adata = basic_qc_plot(axes, adata, ed_low)
        plt.tight_layout()

        if save:
            plt.savefig(save+'/%s_basic_qc.png'%sample, dpi=150, transparent=False, facecolor='white')
            plt.close()
        else:
            plt.show()

    else:

        if filter == "CellBender":

            # Read CellBender results from outside
            cb_out_path = 'cellbender_out/%s/cb_out.h5'%sample
            if not os.path.exists(cb_out_path):
                raise HaltException('No CellBender results, pipeline aborted.')

            with tables.open_file(cb_out_path) as f:
                cb_barcode_indices = f.root.background_removed.barcode_indices_for_latents.read()
                cb_barcodes = adata.obs.index[cb_barcode_indices]
                #cb_barcodes = [i.decode("utf-8") for i in cb_barcodes]
                cb_cell_prob = f.root.background_removed.latent_cell_probability.read()
                f.close()

            adata.obs['cb_prob'] = np.nan
            adata.obs.loc[cb_barcodes, 'cb_prob'] = cb_cell_prob
            adata.obs['cb_class'] = ['Cell' if i > 0.5 else 'Empty droplet' for i in adata.obs['cb_prob']]

            adata_ = sc.read_10x_h5(cb_out_path, genome='background_removed')
            adata_.var = adata.var
            adata_.var.columns = [i+'_original' if '_counts' in i else i for i in adata_.var.columns]
            adata_.obs = adata.obs
            adata_.obs.columns = [i+'_original' if '_counts' in i else i for i in adata_.obs.columns]
            adata_.uns = adata.uns

            sc.pp.calculate_qc_metrics(adata_, qc_vars=['mt'], percent_top=None, log1p=True, inplace=True)
            adata = adata_
        
            # Filtering droplets
            adata = adata[adata.obs.query('total_counts_original > 0').index, adata.var.query('n_cells_by_counts_original > 0').index]
            adata_filtered = adata[adata.obs.query('cb_class == "Cell"').query('pct_counts_mt < 25').index, ]

        elif filter == "EmptyDrops":

            # Filtering droplets
            adata = adata[adata.obs.query('total_counts > 0').index, adata.var.query('n_cells_by_counts > 0').index]
            adata_filtered = adata[adata.obs.query('ed_class == "Cell"').query('pct_counts_mt < 25').index, ]
        
        print('Unfiltered Data dimension: %s droplets x %s genes'%adata.shape)
        print('Pre-filtered Data dimension: %s droplets x %s genes'%adata_filtered.shape)

        # Doublet detection
        print('\nRunning Scrublet...')
        sc.external.pp.scrublet(adata_filtered)
        print('Finished.\n')

        adata_filtered = preprocessing(adata_filtered)
        fig_meta, fig_genes = advanced_qc_plot(adata_filtered, marker_genes, filter=filter)

        if save:
            fig_meta.savefig(save+'/%s_advanced_qc_meta.png'%sample, dpi=150, transparent=False, facecolor='white')
            plt.close(fig_meta)
            fig_genes.savefig(save+'/%s_advanced_qc_marker_expression.png'%sample, dpi=150, transparent=False, facecolor='white')
            plt.close(fig_genes)
        else:
            plt.show()

        # Write output
        print('Writing h5ad output...')
        if not os.path.exists('./processed_anndata'):
            os.mkdir('./processed_anndata')
        adata.write('processed_anndata/%s_unfiltered.h5ad'%sample)
        adata_filtered.write('processed_anndata/%s_filtered_tmp.h5ad'%sample)
        print('Finished.')

def run_pipeline_two(sample, clusters_to_remove=[], ck_use_ref=True, clusters_tumor=[], ck_n_cores=32, save=None, skip_ck=False):

    if os.path.exists('processed_anndata/%s_unfiltered.h5ad'%sample):
        print('Reading pre-filtered data of %s'%sample)
        adata = sc.read_h5ad('processed_anndata/%s_unfiltered.h5ad'%sample)
        adata_filtered = sc.read_h5ad('processed_anndata/%s_filtered_tmp.h5ad'%sample)
        print('Finished.\n')
    else:
        raise HaltException('Pre-filtered AnnData for sample %s do not exist!'%sample)

    filtered_barcodes = adata_filtered.obs.query('leiden not in @clusters_to_remove').index.tolist()
    adata_filtered = adata[filtered_barcodes, ]
    
    # Doublet detection
    print('Running Scrublet...')
    sc.external.pp.scrublet(adata_filtered)
    print('Removing %s doublets...'%adata_filtered.obs['predicted_doublet'].sum())
    adata_filtered = adata_filtered[~adata_filtered.obs['predicted_doublet'], ]
    print('Finished.\n')

    print('Filtered Data dimension: %s droplets x %s genes'%adata_filtered.shape)

    # Print mean stats
    print('\nAverage UMI count:')
    print(adata_filtered.obs['total_counts'].sum()/len(adata_filtered.obs))
    print('Average gene count:')
    print(adata_filtered.obs['n_genes_by_counts'].sum()/len(adata_filtered.obs))

    if not skip_ck:

        # Aneuploidy estimation
        print('\nRunning CopyKat...')
        sample = adata_filtered.uns['sample']
        df_adata_filtered = pd.DataFrame(adata_filtered.X.T.toarray(), index=adata_filtered.var.index, columns=adata_filtered.obs.index)
        df_adata_filtered.index = [i.split('.')[0] for i in df_adata_filtered.index]

        if ck_use_ref:

            normal_cell_barcodes = adata_filtered.obs.query('leiden not in @clusters_tumor').index.tolist()
            normal_cell_barcodes = [i.replace('-', '.') for i in normal_cell_barcodes]

            ck_out = copykat.copykat(df_adata_filtered, id_type="Ensemble", ngene_chr=5, \
                                   win_size=25, KS_cut=0.1, sam_name=sample, distance="euclidean", \
                                   norm_cell_names=normal_cell_barcodes, n_cores=ck_n_cores)
        else:
            ck_out = copykat.copykat(df_adata_filtered, id_type="Ensemble", ngene_chr=5, \
                                   win_size=25, KS_cut=0.1, sam_name=sample, distance="euclidean", \
                                   n_cores=ck_n_cores)
        if not os.path.exists('./copykat_out'):
            os.mkdir('./copykat_out')
        if not os.path.exists('./copykat_out/%s'%sample):
            os.mkdir('./copykat_out/%s'%sample)
        subprocess.call('mv *%s_copykat* ./copykat_out/%s/'%(sample, sample), shell=True)
        
        adata_filtered.obs['ploidy'] = np.nan
        adata_filtered.obs.loc[[i.replace('.', '-') for i in list(ck_out[0][:int(len(ck_out[0])/2)])], 'ploidy'] = list(ck_out[0][int(len(ck_out[0])/2):])
        print('Finished.\n')

    else:

        print('\nSkipping CNV inference.\n')
        adata_filtered.obs['ploidy'] = 'Putatively diploid'

    # Preprocessing
    adata_filtered = preprocessing(adata_filtered)

    with rc_context({'figure.figsize': (4, 4)}):
        fig_meta = sc.pl.embedding(adata_filtered, basis='X_umap', color=['leiden', 'log1p_total_counts', 'pct_counts_mt', 'doublet_score', 'ed_fdr', 'ploidy'], \
                        add_outline=False, legend_loc='right margin',
                        legend_fontsize=8, legend_fontoutline=2, frameon=False,
                        title=['Leiden cluster', 'Total UMI count', 'Mitochondrial UMI rate', 'Doublet score', 'EmptyDrops cell FDR', 'Ploidy'], ncols=3, return_fig=True)
    
    if save:
        fig_meta.savefig(save+'/%s_advanced_qc_meta_with_cnv.png'%sample, dpi=150, transparent=False, facecolor='white')
        plt.close(fig_meta)
    else:
        fig_meta.show()

    # Write output
    print('Writing h5ad output...')
    adata_filtered.write('processed_anndata/%s_filtered.h5ad'%sample)
    print('Finished.')
