import sys
import os
sys.path.append(os.path.abspath('./OV_utility_functions'))
from basic_imports import *
#from VDJseq_analysis_pipeline import *
#from scRNAseq_integrative_analysis_pipeline import *

from itertools import product, combinations
#from more_itertools import unique_everseen

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

import scanpy as sc
#from statannotations.Annotator import Annotator

#import gseapy
#from gseapy.plot import gseaplot

#from upsetplot import UpSet, from_contents

#import anndata2ri
#import rpy2.rinterface_lib.callbacks
#import rpy2.robjects as ro

sc.settings._vector_friendly = True
sc.settings._frameon = False
sc.settings._transparent = True

#def plot_circos():

# Visualize integration quality

def plot_diff_lollipop(x1, x2, label1, label2, xlabel, yticklabels):

    y_range = range(1, len(x1)+1)
    plt.figure(figsize=(6, 7))
    plt.hlines(y=y_range, xmin=x1, xmax=x2, color='grey', alpha=0.6, zorder=-10)
    plt.scatter(x1, y_range, color='skyblue', alpha=1, label=label1)
    plt.scatter(x2, y_range, color='darksalmon', alpha=1 , label=label2)
    plt.legend()
     
    plt.yticks(y_range, yticklabels)
    #plt.title("Comparison of the value 1 and the value 2", loc='left')
    plt.xlabel(xlabel)
    plt.ylabel('')
    sns.despine(right=True, top=True)

def plot_silhouette_avg(adata, df_sil, df_sil_to_contrast, cluster, group, label1, label2, group_order=None):

    # Reformatting silhouette results
    groups = adata.obs[group].unique()
    indices = []
    for i in groups:
        indices.extend(adata.obs[adata.obs[group]==i].index)
    df_sil_ = df_sil.copy().fillna(0)
    df_sil_.index = indices[:len(df_sil_)]
    df_sil_['cluster'] = adata.obs.loc[df_sil_.index, cluster].tolist()

    df_sil_to_contrast_ = df_sil_to_contrast.copy().fillna(0)
    df_sil_to_contrast_['cluster'] = df_sil_['cluster'].tolist()

    df_sil_ = df_sil_.groupby(['cluster', 'group']).mean().reset_index()
    df_sil_['condition'] = label1
    df_sil_to_contrast_ = df_sil_to_contrast_.groupby(['cluster', 'group']).mean().reset_index()
    df_sil_to_contrast_['condition'] = label2
    df_sil_ = pd.concat([df_sil_, df_sil_to_contrast_], ignore_index=True)

    if not group_order:
        group_order = sorted(df_sil_['group'].unique())

    plt.figure(figsize=(12, 6))
    ax = sns.swarmplot(x='group', y='silhouette_score', hue='condition', data=df_sil_, dodge=True, s=3, order=group_order, hue_order=[label1, label2])
    plt.xticks(rotation=45, ha='right')
    sns.despine(left=True, bottom=True, right=True, top=True)
    plt.xlabel('')
    plt.ylabel('Average silhouette coefficient')
    legend = ax.legend([label1, label2], title=None, loc=3)
    legend.get_title().set_size(12)
    for i in range(len(legend.legendHandles)):
        legend.legendHandles[i]._sizes = [50]

    n = -0.15
    for i in group_order:
        x = df_sil_.query('group==@i and condition==@label1').sort_values('cluster')['silhouette_score']
        y = df_sil_.query('group==@i and condition==@label2').sort_values('cluster')['silhouette_score']
        p = st.wilcoxon(x, y)[1]
        if p <= 0.001:
            star = '***'
        elif p <= 0.01:
            star = ' **'
        elif p <= 0.05:
            star = '  *'
        else:
            star = 'ns'
        if x.mean() < y.mean():
            color='k'
        else:
            color='k'
        ax.text(n, df_sil_['silhouette_score'].max()+0.05, star, color=color)
        for j in range(len(x)):
            plt.plot([n-0.05,n+0.35], [x.tolist()[j], y.tolist()[j]], c='k', linewidth=0.3)
        #plt.plot([n-0.15, n+0.1], [y.mean(), y.mean()], c='r', linewidth=2, zorder=10)
        #plt.plot([n+0.2, n+0.45], [x.mean(), x.mean()], c='r', linewidth=2, zorder=10)
        n += 1

def plot_silhouette(adata, df_sil, cluster, group, cluster_order=None, group_subset=None, df_sil_to_contrast=pd.DataFrame({}), label1=None, label2=None):

    # Reformatting silhouette results
    groups = adata.obs[group].unique()
    indices = []
    for i in groups:
        indices.extend(adata.obs[adata.obs[group]==i].index)
    df_sil_ = df_sil.copy()
    df_sil_.index = indices[:len(df_sil_)]
    df_sil_['cluster'] = adata.obs.loc[df_sil_.index, cluster]

    if len(df_sil_to_contrast) > 0:
        df_sil_to_contrast_ = df_sil_to_contrast.copy()
        df_sil_to_contrast_['cluster'] = df_sil_['cluster'].tolist()

    if not cluster_order:
        cluster_order = list(df_sil_['cluster'].unique())

    if group_subset:
        df_sil_ = df_sil_.query('group == @group_subset')
        if len(df_sil_to_contrast) > 0:
            df_sil_to_contrast_ = df_sil_to_contrast_.query('group == @group_subset')

    n_clusters = len(df_sil_['cluster'].unique())

    interval = int(len(df_sil_)/50)

    fig, ax = plt.subplots(1, 1, figsize=(7, 7))
    if len(df_sil_to_contrast) > 0:
        xlim = min(df_sil_['silhouette_score'].min(), df_sil_to_contrast_['silhouette_score'].min())-0.02, max(df_sil_['silhouette_score'].max(), df_sil_to_contrast_['silhouette_score'].max())+0.02
        
    else:
        xlim = df_sil_['silhouette_score'].min()-0.02, df_sil_['silhouette_score'].max()+0.02
    ax.set_xlim(xlim)
    ax.set_ylim([0, len(df_sil_) + (n_clusters + 1) * interval])
    y_lower = interval
    for cluster in cluster_order:
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        i = cluster_order.index(cluster)
        ith_cluster_silhouette_values = df_sil_.query('cluster == @cluster')['silhouette_score'].sort_values(ascending=True)
        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i
        color = plt.cm.Spectral_r(float(i) / n_clusters)
        ax.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor='k',
            alpha=1,
        )

        if len(df_sil_to_contrast) > 0:
            ith_cluster_silhouette_values_to_contrast = df_sil_to_contrast_.query('cluster == @cluster')['silhouette_score'].sort_values(ascending=True)
            ax.fill_betweenx(
                np.arange(y_lower, y_upper),
                ith_cluster_silhouette_values,
                ith_cluster_silhouette_values_to_contrast,
                facecolor=color,
                edgecolor='k',
                alpha=0.5,
                linestyle='--'
            )

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(xlim[0]-0.02, y_lower + 0.5 * size_cluster_i, cluster, ha='right', va='center')

        # Compute the new y_lower for next plot
        y_lower = y_upper + interval  # 10 for the 0 samples

    ax.set_xlabel("Silhouette coefficient")

    # The vertical line for average silhouette score of all the values
    if len(df_sil_to_contrast) > 0:
        silhouette_avg = df_sil_['silhouette_score'].mean()
        ax.axvline(x=silhouette_avg, color="red", linestyle="-", lw=1, label=label1)
        silhouette_avg_to_contrast = df_sil_to_contrast_['silhouette_score'].mean()
        ax.axvline(x=silhouette_avg_to_contrast, color="red", linestyle="--", lw=1, alpha=0.5, label=label2)
    else:
        silhouette_avg = df_sil_['silhouette_score'].mean()
        ax.axvline(x=silhouette_avg, color="red", linestyle="--", lw=1, label='Average')

    ax.set_yticks([])  # Clear the yaxis labels / ticks

    sns.despine(ax=ax, right=True, top=True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles[::-1], labels[::-1], title='', loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6}, frameon=False)

# Visualize sample-wise feature association

def linear_percentile(values, perc):
    return (perc/100)*(max(values) - min(values))+min(values)

def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

def hanging_line(point1, point2):
    import numpy as np

    a = (point2[1] - point1[1])/(np.cosh(point2[0]) - np.cosh(point1[0]))
    b = point1[1] - a*np.cosh(point1[0])
    x = np.linspace(point1[0], point2[0], 100)
    y = a*np.cosh(x) + b

    return (x,y)

def fancy_heatmap(fig, ax, df, df_size=None, only_size=False, unit='count', **kwargs):

    df_long = df.T
    df_long.index.name = 'index'
    df_long = df_long.reset_index()
    df_long = pd.melt(df_long, id_vars=['index'], value_vars=df_long.columns[1:])

    x = df_long['index']
    y = df_long['variable']
    z = df_long['value']

    if 'palette' in kwargs:
        palette = kwargs['palette']
        n_colors = len(palette)
    else:
        n_colors = 256 # Use 256 colors for the diverging color palette
        palette = sns.color_palette("Blues", n_colors) 

    if only_size:
        color = [a for b in [[i]*df.shape[1] for i in range(df.shape[0])] for a in b]
        ax, ax_slgd = ax
    else:
        color = z
        ax, ax_slgd, ax_cbar = ax
    
    try:
        df_size_long = df_size.T
        df_size_long = df_size_long.reset_index()
        df_size_long = pd.melt(df_size_long, id_vars=['index'], value_vars=df_size_long.columns[1:])
        size = df_size_long['value']
    except:
        size = z
        
    if 'color_range' in kwargs:
        color_min, color_max = kwargs['color_range']
    else:
        color_min, color_max = min(color), max(color) # Range of values that will be mapped to the palette, i.e. min and max possible correlation

    def value_to_color(val):
        if color_min == color_max:
            return palette[-1]
        else:
            val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            ind = int(val_position * (n_colors - 1)) # target index in the color palette
            return palette[ind]

    if 'size_range' in kwargs:
        size_min, size_max = kwargs['size_range'][0], kwargs['size_range'][1]
    else:
        size_min, size_max = min(size), max(size)

    size_scale = kwargs.get('size_scale', 500)

    def value_to_size(val):
        if size_min == size_max:
            return 1 * size_scale
        else:
            val_position = (val - size_min) * 0.99 / (size_max - size_min) # position of value in the input range, relative to the length of the input range
            val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
            return val_position * size_scale

    if 'x_order' in kwargs: 
        x_names = [t for t in kwargs['x_order']]
    else:
        x_names = [t for t in sorted(set([v for v in x]))]
    x_to_num = {p[1]:p[0] for p in enumerate(x_names)}

    if 'y_order' in kwargs: 
        y_names = [t for t in kwargs['y_order']][::-1]
    else:
        y_names = [t for t in sorted(set([v for v in y]))][::-1]
    y_to_num = {p[1]:p[0] for p in enumerate(y_names)}

    marker = kwargs.get('marker', 's')

    kwargs_pass_on = {k:v for k,v in kwargs.items() if k not in [
         'color', 'palette', 'color_range', 'size', 'size_range', 'size_scale', 'marker', 'x_order', 'y_order', 'xlabel', 'ylabel'
    ]}

    ax.scatter(
        x=[x_to_num[v] for v in x],
        y=[y_to_num[v] for v in y],
        marker=marker,
        s=[value_to_size(v) for v in size], 
        c=[value_to_color(v) for v in color],
        **kwargs_pass_on
        )

    ax.set_xticks([v for k,v in x_to_num.items()])
    ax.set_xticklabels([k for k in x_to_num], rotation=45, horizontalalignment='right')
    ax.set_yticks([v for k,v in y_to_num.items()])
    ax.set_yticklabels([k for k in y_to_num])

    ax.grid(False, 'major')
    ax.grid(False, 'minor')

    ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
    ax.set_ylim([-0.5, max([v for v in y_to_num.values()]) + 0.5])
    #ax.set_facecolor('#F1F1F1')

    ax.set_xlabel(kwargs.get('xlabel', ''))
    ax.set_ylabel(kwargs.get('ylabel', ''))

    # Add size legend on the bottom left region of the plot
    ax_slgd_num = 4
    ax_slgd_xlim_extend = 5
    ax_slgd_width = get_ax_size(fig, ax_slgd)[0]
    ax_slgd_dpi = fig.dpi
    ax_slgd_spacing = 0.3
    ax_slgd_single_data_unit_length_in_point = (ax_slgd_width  / (ax_slgd_num+ax_slgd_xlim_extend+1) * 72./350)
    ax_slgd_sizes_in_point = [linear_percentile([value_to_size(v) for v in size], perc) for perc in np.arange(5, 100+1, 100/(ax_slgd_num))]
    ax_slgd_lengths_in_point = [np.sqrt(i) for i in ax_slgd_sizes_in_point]
    ax_slgd_lengths_in_data_unit = [i/ax_slgd_single_data_unit_length_in_point for i in ax_slgd_lengths_in_point]
    def get_ax_slgd_center(n):
        if n == 0:
            return 0
        else:
            return get_ax_slgd_center(n-1) + ax_slgd_lengths_in_data_unit[n-1]/2 + ax_slgd_spacing + ax_slgd_lengths_in_data_unit[n]/2
    ax_slgd_centers = [get_ax_slgd_center(i) for i in range(ax_slgd_num)]

    if unit=='proportion':
        ax_slgd_labels = ['%.0f'%linear_percentile(size, perc) for perc in np.arange(5, 100+1, 100/(ax_slgd_num))]
        ax_slgd_title = 'Proportion (%)'
    elif unit=='count':
        ax_slgd_labels = ['%.1e'%linear_percentile(size, perc) for perc in np.arange(5, 100+1, 100/(ax_slgd_num))]
        ax_slgd_title = 'Count'
    elif unit=='log_count':
        ax_slgd_labels = ['%.1f'%linear_percentile(size, perc) for perc in np.arange(5, 100+1, 100/(ax_slgd_num))]
        ax_slgd_title = 'Count (log10)'
    else:
        print('Not supported')

    ax_slgd.axis('off')
    ax_slgd.set_xlim(-1, ax_slgd_num+ax_slgd_xlim_extend)
    ax_slgd.set_ylim(0, 1)
    ax_slgd.scatter(
            x=ax_slgd_centers,
            y=[0.5]*ax_slgd_num,
            marker=marker,
            s=ax_slgd_sizes_in_point, 
            c='k'
            )
    [ax_slgd.text(ax_slgd_centers[i], 0.5-max(ax_slgd_lengths_in_data_unit)/2/(ax_slgd_num+5+1)-0.05, ax_slgd_labels[i], ha='center', va='top', fontsize=10) for i in range(ax_slgd_num)]
    ax_slgd.text(min(ax_slgd_centers)+(max(ax_slgd_centers)-min(ax_slgd_centers))/2, 0.5+max(ax_slgd_lengths_in_data_unit)/2/(ax_slgd_num+5+1)+0.05, ax_slgd_title, ha='center', va='bottom', fontsize=12)

    if not only_size:
        # Add color legend on the right side of the plot

        col_x = [0]*len(palette) # Fixed x coordinate for the bars
        bar_y=np.linspace(color_min, color_max, n_colors) # y coordinates for each of the n_colors bars

        bar_height = bar_y[1] - bar_y[0]
        ax_cbar.barh(
            y=bar_y,
            width=[5]*len(palette), # Make bars 5 units wide
            left=col_x, # Make bars start at 0
            height=bar_height,
            color=palette,
            linewidth=0
        )
        
        ax_cbar.set_xlim(1, 2) # Bars are going from 0 to 5, so lets crop the plot somewhere in the middle
        #ax_cbar.grid(False, 'major') # Hide grid
        #ax_cbar.grid(False, 'minor') # Hide grid
        [ax_cbar.spines[i].set_visible(False) for i in ['left', 'right', 'top', 'bottom']]
        ax_cbar.set_facecolor('white') # Make background white
        ax_cbar.set_xticks([]) # Remove horizontal ticks
        ax_cbar.set_yticks(np.linspace(min(bar_y), max(bar_y), 3)) # Show vertical ticks for min, middle and max
        ax_cbar.set_ylabel(ax_slgd_title)
        ax_cbar.yaxis.set_label_position("right")
        ax_cbar.yaxis.tick_right() # Show vertical ticks on the right

def fancy_barplot(ax, df, legend=True, bar_cmap=palettable.cartocolors.qualitative.Prism_10.mpl_colors, bar_width=0.6, bar_linkage=True, bar_link_alpha=0.5):

    colors = [bar_cmap[i] for i in range(len(df))][::-1]
    df = df.loc[df.index[::-1], ]

    [ax.bar(\
         x=df.columns, height=df.iloc[i, :], \
         bottom = df.iloc[:i, :].sum(axis=0).tolist(), \
         label = df.index[i], color = colors[i], width=bar_width \
        ) for i in range(len(df.index))]

    if legend:
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
        handles, labels = ax.get_legend_handles_labels()
        legend = ax.legend(handles[::-1], labels[::-1], title='', loc='center left', bbox_to_anchor=(1, 0.5), prop={'size': 6}, frameon=False)
        #legend.set_title(y.title())
        legend.get_title().set_size(7)

    ax.set_xticklabels(df.columns, rotation=45, ha='right', fontsize=8)
    sns.despine(ax=ax, top=True, right=True)

    if bar_linkage:
        nrows, ncols = df.shape
        for nrow in range(nrows):
            for ncol in range(ncols-1):
                if nrow == 0:
                    x1, y1 = hanging_line([ncol+bar_width/2, 0], \
                                          [ncol+1-(bar_width/2), 0])
                else:
                    x1, y1 = hanging_line([ncol+bar_width/2, df.iloc[:nrow, ncol].sum()], \
                                          [ncol+1-(bar_width/2), df.iloc[:nrow, ncol+1].sum()])
                x2, y2 = hanging_line([ncol+bar_width/2, df.iloc[:nrow+1, ncol].sum()], \
                                      [ncol+1-(bar_width/2), df.iloc[:nrow+1, ncol+1].sum()])
                ax.fill_between(x1, y1, y2, color=colors[nrow], alpha=bar_link_alpha, edgecolor=None)

def plot_composition_simple(adata, x, y, x_order=None, y_order=None, norm=True, style='bar', colors=None, cmap=palettable.cartocolors.qualitative.Prism_10.get_mpl_colormap(), \
                            bar_figsize=(8, 6), bar_width=0.6, bar_linkage=True, bar_link_alpha=0.5, \
                            heat_figsize=(5, 5), only_size=False, heat_size_scale=50, heat_marker='.', \
                            return_df=False, melt=True, return_df_annot=None):

    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[y].unique()))

    if not x_order:
        x_order = df[x].unique()

    if not y_order:
        y_order = df_frac.index

    for i in x_order:
        df_frac[i] = df[df[x]==i].groupby(y).count().iloc[:, 0]

    df_frac = df_frac.loc[y_order, :]

    if norm:
        df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)*100
        unit = 'proportion'
        bar_ylabel = 'Proportion (%)'
    else:
        df_frac = df_frac+1
        unit = 'count'
        bar_ylabel = 'Count'

    if colors:
        cmap = colors
    else:
        cmap = [plt.cm.get_cmap(cmap)(perc/100) for perc in np.arange(0, 100, 100/len(y_order))]

    if style == 'bar':

        fig, ax = plt.subplots(1, 1, figsize=(bar_figsize))

        fancy_barplot(ax=ax, df=df_frac, bar_cmap=cmap, bar_width=bar_width, bar_linkage=bar_linkage, bar_link_alpha=bar_link_alpha)
        ax.set_ylabel(bar_ylabel)

    elif style == 'heatmap':

        fig = plt.figure(figsize=heat_figsize)
        if only_size:
            plot_grid = plt.GridSpec(5, 5, hspace=0.2, wspace=0.2)
            ax = (plt.subplot(plot_grid[:4,1:-1]), plt.subplot(plot_grid[4:,:1]))
        else:
            plot_grid = plt.GridSpec(4, 16, hspace=0.2, wspace=0.2)
            ax = (plt.subplot(plot_grid[:3,3:-1]), plt.subplot(plot_grid[3,:3]), plt.subplot(plot_grid[:3,-1]))

        fancy_heatmap(fig=fig, ax=ax, df=df_frac, only_size=only_size, unit=unit, 
                      x_order=x_order, y_order=y_order, 
                      size_scale=heat_size_scale, palette=cmap, marker=heat_marker)

    if return_df:
        if melt:
            df_frac_long = df_frac.T
            df_frac_long = df_frac_long.reset_index()
            df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
            df_frac_long = df_frac_long.rename(columns={'index':x, 'variable':y, 'value':'fraction'})
            if return_df_annot:
                for i in return_df_annot:
                    df_frac_long[i] = adata.obs.drop_duplicates(x).set_index(x).loc[df_frac_long[x], i].tolist()
            return df_frac_long
        else:
            df_frac = df_frac.T
            if return_df_annot:
                for i in return_df_annot:
                    df_frac[i] = adata.obs.drop_duplicates(x).set_index(x).loc[df_frac.index, i].tolist()
            return df_frac

def plot_composition_simple2(adata, x, y, x_order=None, y_order=None, norm=True, style='bar', colors=None, cmap=palettable.cartocolors.qualitative.Prism_10.get_mpl_colormap(), \
                            bar_figsize=(8, 6), bar_width=0.6, bar_linkage=True, bar_link_alpha=0.5, \
                            heat_figsize=(5, 5), only_size=False, heat_size_scale=50, heat_marker='.', \
                            return_df=False, melt=True, return_df_annot=None):

    # Parse data
    df = adata.obs
    df_frac_tmp = adata.obsm[y]
    df_frac = pd.DataFrame({}, index=df_frac_tmp.columns)

    if not x_order:
        x_order = df[x].unique()

    if not y_order:
        y_order = df_frac_tmp.columns

    for i in x_order:
        df_frac[i] = df_frac_tmp.loc[df[x]==i,:].mean(axis=0)
    df_frac = df_frac.loc[y_order, :]

    if norm:
        df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)*100
        unit = 'proportion'
        bar_ylabel = 'Proportion (%)'
    else:
        df_frac = df_frac
        unit = 'proportion'
        bar_ylabel = 'proportion'

    if colors:
        cmap = colors
    else:
        cmap = [plt.cm.get_cmap(cmap)(perc/100) for perc in np.arange(0, 100, 100/len(y_order))]

    if style == 'bar':
        fig, ax = plt.subplots(1, 1, figsize=(bar_figsize))
        fancy_barplot(ax=ax, df=df_frac, bar_cmap=cmap, bar_width=bar_width, bar_linkage=bar_linkage, bar_link_alpha=bar_link_alpha)
        ax.set_ylabel(bar_ylabel)

    elif style == 'heatmap':
        fig = plt.figure(figsize=heat_figsize)
        if only_size:
            plot_grid = plt.GridSpec(5, 5, hspace=0.2, wspace=0.2)
            ax = (plt.subplot(plot_grid[:4,1:-1]), plt.subplot(plot_grid[4:,:1]))
        else:
            plot_grid = plt.GridSpec(4, 16, hspace=0.2, wspace=0.2)
            ax = (plt.subplot(plot_grid[:3,3:-1]), plt.subplot(plot_grid[3,:3]), plt.subplot(plot_grid[:3,-1]))
        fancy_heatmap(fig=fig, ax=ax, df=df_frac, only_size=only_size, unit=unit, 
                      x_order=x_order, y_order=y_order, 
                      size_scale=heat_size_scale, palette=cmap, marker=heat_marker)

    if return_df:
        if melt:
            df_frac_long = df_frac.T
            df_frac_long = df_frac_long.reset_index()
            df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
            df_frac_long = df_frac_long.rename(columns={'index':x, 'variable':y, 'value':'fraction'})
            if return_df_annot:
                for i in return_df_annot:
                    df_frac_long[i] = adata.obs.drop_duplicates(x).set_index(x).loc[df_frac_long[x], i].tolist()
            return df_frac_long
        else:
            df_frac = df_frac.T
            if return_df_annot:
                for i in return_df_annot:
                    df_frac[i] = adata.obs.drop_duplicates(x).set_index(x).loc[df_frac.index, i].tolist()
            return df_frac

def plot_composition_complex(adata, x, y, x_order=None, y_order=None, x_annot=None, xticklabels=True, \
                             figsize=(10, 8), colors=None, palette=palettable.cartocolors.qualitative.Prism_10.get_mpl_colormap(), \
                             heat_norm=True, bar_x_norm=False, bar_y_norm=True, \
                             bar_width=0.6, bar_x_linkage=True, bar_x_link_alpha=0.5, \
                             heat_size_scale=50, heat_marker='.'):

    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[y].unique()))

    if not x_order:
        x_order = list(df[x].unique())

    if not y_order:
        y_order = df_frac.index

    for i in x_order:
        df_frac[i] = df[df[x]==i].groupby(y).count().iloc[:, 0]

    df_frac = df_frac.loc[y_order, :]
    df_frac_heat = df_frac.copy()
    df_frac_bar_x = df_frac.copy()
    df_frac_bar_y = df_frac.copy()

    if heat_norm:
        df_frac_heat = df_frac_heat.divide(df_frac_heat.sum(axis=0), axis=1)*100
        heat_unit = 'proportion'
    else:
        df_frac_heat = df_frac_heat
        heat_unit = 'count'

    if bar_x_norm:
        df_frac_bar_x = df_frac_bar_x.divide(df_frac_bar_x.sum(axis=0), axis=1)*100
        bar_x_ylabel = 'Proportion (%)'
    else:
        df_frac_bar_x = df_frac_bar_x
        bar_x_ylabel = 'Count'

    if bar_y_norm:
        vec_frac_bar_y = df_frac_bar_y.sum(axis=1)/df_frac_bar_y.sum(axis=0).sum()*100
        bar_y_xlabel = 'Proportion (%)'
    else:
        vec_frac_bar_y = np.log10(df_frac_bar_y.sum(axis=1)+1)
        bar_y_xlabel = 'Count (log10)'

    if not colors:
        palette = [plt.cm.get_cmap(palette)(perc/100) for perc in np.arange(0, 100, 100/len(y_order))]
    else:
        palette = colors

    # Create Fig and gridspec
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(nrows=4, ncols=4, hspace=0.05, wspace=0.05)

    # Define the axes
    ax_heatmap_main = fig.add_subplot(grid[1:3, 1:3])
    ax_heatmap_main_legend = fig.add_subplot(grid[3, 0])
    if x_annot:
        ax_heatmap_x_annot = fig.add_subplot(grid[3, 1:3])
    ax_heatmap_y_annot = fig.add_subplot(grid[1:3, 0])
    ax_bar_x = fig.add_subplot(grid[0, 1:3])
    ax_bar_y = fig.add_subplot(grid[1:3, 3])

    # Configure the axes
    ax_heatmap_main.set_zorder(10)
    ax_heatmap_main.tick_params(axis='y', which='major', pad=15)
    #ax_heatmap_main.add_patch(Rectangle((-0.5, -0.5), 8, 8,
    #                  alpha=1, fill=None, linewidth=2))

    #ax_heatmap_y_annot.margins(y=0.01)
    ax_heatmap_y_annot.set_ylim(-0.5, len(df_frac_bar_y.index)-0.5)
    ax_heatmap_y_annot.axis('off')

    #ax_bar_x.margins(x=0.01)
    ax_bar_x.set_ylabel(bar_x_ylabel)
    ax_bar_x.set_xlim(-0.5, len(df_frac_bar_x.columns)-0.5)
    ax_bar_x.get_xaxis().set_visible(False)
    ax_bar_x.locator_params(nbins=5, axis='y')
    sns.despine(ax=ax_bar_x, left=False, right=True, top=True, bottom=True)

    ax_bar_y.get_yaxis().set_visible(False)
    sns.despine(ax=ax_bar_y, left=True, right=True, top=False, bottom=True)
    ax_bar_y.locator_params(nbins=5, axis='x')
    ax_bar_y.set_xlabel(bar_y_xlabel)
    ax_bar_y.xaxis.set_label_position("top")
    ax_bar_y.xaxis.tick_top()

    # Draw the main heatmap
    fancy_heatmap(fig=fig, ax=(ax_heatmap_main, ax_heatmap_main_legend), df=df_frac_heat, only_size=True, unit=heat_unit, \
                  x_order=x_order, y_order=y_order, \
                  size_scale=heat_size_scale, palette=palette, marker=heat_marker)
    if not xticklabels:
        ax_heatmap_main.set_xticks([])

    color = list(range(len(y_order)))
    n_colors = len(palette)
    color_min, color_max = min(color), max(color) # Range of values that will be mapped to the palette, i.e. min and max possible correlation
    def value_to_color(val):
        val_position = float((val - color_min)) / (color_max - color_min) # position of value in the input range, relative to the length of the input range
        val_position = min(max(val_position, 0), 1) # bound the position betwen 0 and 1
        ind = int(val_position * (n_colors - 1)) # target index in the color palette
        return palette[ind]
    palette = [value_to_color(v) for v in color]

    # Draw the left-sided annotation heatmap
    ax_heatmap_y_annot.barh(y=df_frac_bar_y.index, width=[1]*len(df_frac_bar_y.index), \
                edgecolor=[palette[i] for i in range(len(vec_frac_bar_y))][::-1], alpha=1, height=0.9, fill=False, linewidth=1, linestyle='--')
    ax_heatmap_y_annot.scatter(x=[0.1]*len(df_frac_bar_y.index), y=list(range(len(df_frac_bar_y.index))), \
                color=palette[::-1], alpha=1)

    # Draw the bottom annotation heatmap
    if x_annot:
        df_x_annot = adata.obs.drop_duplicates(x).set_index(x).loc[x_order, x_annot]
        for i in x_annot:
            if isinstance(df_x_annot[i][0], str):
                df_x_annot[i] = df_x_annot[i].map(dict(zip(df_x_annot[i].unique(), range(len(df_x_annot[i].unique()))))).tolist()
        df_x_annot = df_x_annot.T
        df_x_annot = st.zscore(df_x_annot, axis=1)
        sns.heatmap(df_x_annot, ax=ax_heatmap_x_annot, linewidth=0.5, cbar=False, cmap=plt.cm.Spectral_r)
        ax_heatmap_main.get_xaxis().set_visible(False)
        ax_heatmap_x_annot.set_xlabel('')
        ax_heatmap_x_annot.yaxis.set_label_position("right")
        ax_heatmap_x_annot.yaxis.tick_right()
        if xticklabels:
            ax_heatmap_x_annot.set_xticks(list(range(len(x_order))), x_order, rotation=45, ha='right', va='top')
        else:
            ax_heatmap_x_annot.set_xticks([])
        ax_heatmap_x_annot.set_yticks([i+0.5 for i in list(range(len(x_annot)))], x_annot, rotation=0, ha='left', va='center')

    # Draw the top stacked bar plot
    fancy_barplot(ax=ax_bar_x, df=df_frac_bar_x, legend=False, \
                  bar_cmap=palette, bar_width=bar_width, bar_linkage=bar_x_linkage, bar_link_alpha=bar_x_link_alpha)

    # Draw the right-sided unstacked bar plot
    sns.barplot(ax=ax_bar_y, x=vec_frac_bar_y, y=df_frac_bar_y.index, \
                palette=palette, orient='h')


def plot_diff_abundance_all(adata, cell_type_col, sample_col,
                            variable, cell_type_order=[], variable_order=[], variable_contrasts=None, 
                            test='t-test_ind', covariate=None, covariate_subset=None, sample_dict=None,
                            colors=plt.cm.Set2, figsize=None, print_pval=False, box_pairs=None):
    # Parse data
    df = adata.obs.copy()
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)

    df_frac_long = df_frac.T
    df_frac_long = df_frac_long.reset_index()
    df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
    df_frac_long = df_frac_long.rename(columns={'index':'Sample', 'variable':'Cell type', 'value':'fraction'})

    if not sample_dict:

        df_frac_long[variable] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac_long['Sample'], variable].tolist()

        if covariate:
            df_frac_long[covariate] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac_long['Sample'], covariate].tolist()
            df_frac_long = df_frac_long[df_frac_long[covariate]==covariate_subset]

    else:

        df_frac_long = df_frac_long[df_frac_long['Sample'].isin([a for b in list(sample_dict.values()) for a in b])]
        df_frac_long[variable] = np.nan
        for i in list(sample_dict.keys()):
            df_frac_long.loc[df_frac_long['Sample'].isin(sample_dict[i]), variable] = i

    if len(cell_type_order) == 0:

        cell_type_order = sorted(df_frac_long['Cell type'].unique())

    if len(variable_order) == 0:

        variable_order = sorted(df_frac_long[variable].unique())

    df_frac_long = df_frac_long[df_frac_long['Cell type'].isin(cell_type_order)]
    df_frac_long = df_frac_long[df_frac_long[variable].isin(variable_order)]

    if not figsize:

        figsize = (len(df_frac_long['Cell type'].unique()), 6)

    if not isinstance(colors, list):

        colors = [colors(i) for i in range(len(variable_order))]

    fig, ax = plt.subplots(1, 1, figsize=(figsize))
    sns.swarmplot(ax=ax, data=df_frac_long, x='Cell type', order=cell_type_order, 
                  y='fraction', hue=variable, hue_order=variable_order, dodge=True, palette=colors)
    sns.boxplot(ax=ax, data=df_frac_long, x='Cell type', order=cell_type_order, 
                y='fraction', hue=variable, hue_order=variable_order, showfliers=False, notch=False,
                whiskerprops={'ls':'--', 'lw':1, 'color':'k'}, 
                capprops={'ls':'--', 'lw':1, 'color':'k'}, 
                boxprops={'ls':'--', 'lw':1, 'edgecolor':'k', 'facecolor':'white'}, 
                medianprops={'ls':'--', 'lw':1, 'color':'k'})
        
    if variable_contrasts or len(variable_order) <= 3:
    
        if print_pval:

            text_format='simple'

        else:

            text_format='star'

        variable_contrasts = combinations(variable_order, 2) if not variable_contrasts else variable_contrasts
        
        box_pairs = [a for b in 
             [[((cell, variable_contrast[0]), 
             (cell, variable_contrast[1])) 
             for cell in df_frac_long['Cell type'].unique()]
             for variable_contrast in variable_contrasts]
             for a in b]
        
        annotator = Annotator(ax=ax, pairs=box_pairs, data=df_frac_long, x='Cell type', order=cell_type_order, 
                              y='fraction', hue=variable, hue_order=variable_order)

        annotator.configure(test=test, comparisons_correction=None, text_format=text_format, loc='inside', show_test_name=False, verbose=0)
        annotator.apply_and_annotate()
    
    else:

        n = -0.15

        val_max = max(df_frac_long['fraction'])
        val_min = min(df_frac_long['fraction'])
        val_range = val_max - val_min

        if not box_pairs:
        
            for i in cell_type_order:
                
                values = [df_frac_long[(df_frac_long['Cell type']==i) & (df_frac_long[variable]==j)]['fraction'].tolist() for j in variable_order]
                
                if test == 'anova':
                    p = st.f_oneway(*values)[1]
                elif test == 'kruskal':
                    p = st.kruskal(*values)[1]
                    
                if p <= 0.001:
                    star = '***'
                elif p <= 0.01:
                    star = '**'
                elif p <= 0.05:
                    star = '*'
                else:
                    star = 'ns'

                if print_pval:
                    ax.text(n+0.15, val_max+val_range*0.02, 'P = %.2g'%p, ha='center')
                else:
                    ax.text(n+0.15, val_max+val_range*0.02, star, ha='center')
                
                n += 1
            
        else:

            if print_pval:

                text_format='simple'

            else:

                text_format='star'

            annotator = Annotator(ax=ax, pairs=box_pairs, data=df_frac_long, x='Cell type', order=cell_type_order, 
                                    y='fraction', hue=variable, hue_order=variable_order)

            annotator.configure(test=test, comparisons_correction=None, text_format=text_format, loc='inside', show_test_name=False, verbose=0)
            annotator.apply_and_annotate()

    lgd = ax.legend(variable_order, title=' '.join(variable.split('_')), frameon=True, prop={'size':11})
    lgd.get_title().set_size(12)
    
    sns.despine(right=True, top=True)
    plt.ylabel('Cell type fraction')
    plt.xlabel('')
    plt.xticks(rotation=45, ha='right')

def plot_diff_abundance_specific(adata, cell_type_col, cell_type, sample_col, 
                                 variable1, variable1_order, variable2, variable2_order, variable1_contrasts=None, 
                                 covariate=None, covariate_subset=None, 
                                 figsize=(5, 5), style='bar', test='t-test_ind', text_format='star', colors=['lightblue', 'darksalmon']):

    # Parse data
    if covariate:
        df = adata.obs.loc[adata.obs[[variable1, variable2, covariate]].dropna().index, ]
    else:
        df = adata.obs.loc[adata.obs[[variable1, variable2]].dropna().index, ]

    # Parse data
    df = adata.obs.copy()
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)

    df_frac_long = df_frac.T
    df_frac_long = df_frac_long.reset_index()
    df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
    df_frac_long = df_frac_long.rename(columns={'index':'Sample', 'variable':'Cell type', 'value':'fraction'})


    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)

    df_frac_long = df_frac.T
    df_frac_long = df_frac_long.reset_index()
    df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
    df_frac_long = df_frac_long.rename(columns={'index':'Sample', 'variable':'Cell type', 'value':'fraction'})

    if isinstance(cell_type, str):

        df_frac_long = df_frac_long[df_frac_long['Cell type'] == cell_type]

    elif len(cell_type) == 1:

        df_frac_long_ = df_frac_long.drop_duplicates('Sample').sort_values('Sample')

        df_frac_long_['fraction'] = df_frac_long[df_frac_long['Cell type'].isin(cell_type[0])].groupby('Sample').sum()['fraction'].values

        df_frac_long = df_frac_long_

    elif len(cell_type) == 2:

        df_frac_long_ratio = df_frac_long.drop_duplicates('Sample').sort_values('Sample')

        df_frac_long_ratio['fraction'] = df_frac_long[df_frac_long['Cell type'].isin(cell_type[0])].groupby('Sample').sum()['fraction'].values / \
                                          df_frac_long[df_frac_long['Cell type'].isin(cell_type[1])].groupby('Sample').sum()['fraction'].values
        df_frac_long_ratio['fraction'] = np.log2(df_frac_long_ratio['fraction'])

        df_frac_long = df_frac_long_ratio

    df_frac_long[variable1] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac_long['Sample'], variable1].tolist()
    df_frac_long[variable2] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac_long['Sample'], variable2].tolist()

    if covariate:
        df_frac_long[covariate] = adata.obs.drop_duplicates(sample_col).set_index(sample_col).loc[df_frac_long['Sample'], covariate].tolist()
        df_frac_long = df_frac_long[df_frac_long[covariate]==covariate_subset]

    if style == 'bar':

        fig, ax = plt.subplots(1, 1, figsize=figsize)
        sns.stripplot(ax=ax, data=df_frac_long, x=variable2, order=variable2_order, 
                      y='fraction', hue=variable1, hue_order=variable1_order, dodge=True, palette = colors)
        sns.boxplot(ax=ax, data=df_frac_long, x=variable2, order=variable2_order, 
                    y='fraction', hue=variable1, hue_order=variable1_order, showfliers=False, notch=False,
                    whiskerprops={'ls':'--', 'lw':1, 'color':'k'}, 
                    capprops={'ls':'--', 'lw':1, 'color':'k'}, 
                    boxprops={'ls':'--', 'lw':1, 'edgecolor':'k', 'facecolor':'white'}, 
                    medianprops={'ls':'--', 'lw':1, 'color':'k'})

        if not variable1_contrasts:
            if len(variable1_order) == 2:
                variable1_contrasts = variable1_order
            else:
                raise KeyError('A variable contrast list has to be input.')
        
        box_pairs = [((i, variable1_contrasts[0]), (i, variable1_contrasts[1])) for i in df_frac_long[variable2].unique()]

        annotator = Annotator(ax=ax, pairs=box_pairs, data=df_frac_long, x=variable2, order=variable2_order, 
                                y='fraction', hue=variable1, hue_order=variable1_order)

        annotator.configure(test=test, comparisons_correction=None, text_format=text_format, loc='inside', show_test_name=False, verbose=0)
        annotator.apply_and_annotate()
        
        lgd = ax.legend(variable1_order, title=' '.join(variable1.split('_')), frameon=True, prop={'size':11})
        lgd.get_title().set_size(12)
        
        sns.despine(right=True, top=True)
        #plt.ylabel('%s fraction'%cell_type)
        plt.xlabel('')
        plt.xticks(rotation=45, ha='right')

    elif style == 'reg':
        
        tmp = [st.spearmanr(df_frac_long[df_frac_long[variable2]==i][variable1], 
                                  df_frac_long[df_frac_long[variable2]==i]['fraction'], nan_policy='omit') for i in variable2_order]
        
        corrs, ps = [i[0] for i in tmp], [i[1] for i in tmp]

        sns.lmplot(data=df_frac_long, x=variable1, y='fraction', 
                   hue=variable2, hue_order=variable2_order, 
                   legend=False, legend_out=True, size=6,
                   palette=colors, robust=True)

        sns.despine(right=True, top=True)
        #plt.ylabel('%s fraction'%cell_type)
        plt.xlabel(' '.join(variable1.split('_')))
        #plt.xticks(rotation=45, ha='right')
        
        handles, labels = plt.gca().get_legend_handles_labels()
        lgd = plt.legend(handles, 
                         [i + '\ncorr=%.2f, p=%.1e'%(corrs[labels.index(i)], ps[labels.index(i)]) for i in labels], 
                         title=' '.join(variable2.split('_')), frameon=True, prop={'size':11}, 
                         bbox_to_anchor=(1, 1))
        lgd.get_title().set_size(12)
            
def plot_diff_abundance_paired(adata, cell_type_col, sample_col, variable, variable_order,
                               sample_dict, sample_pairing_dict, 
                               cell_type_order=[], test='Wilcoxon', colors=['lightblue', 'darksalmon'], figsize=(12, 6), print_pval=False):

    # Parse data
    df = adata.obs.copy()

    if len(cell_type_order) == 0:
        cell_type_order = sorted(df[cell_type_col].unique())
        
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)

    df_frac_long = df_frac.T
    df_frac_long = df_frac_long.reset_index()
    df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
    df_frac_long = df_frac_long.rename(columns={'index':'Sample', 'variable':'Cell type', 'value':'fraction'})

    df_frac_long = df_frac_long[df_frac_long['Sample'].isin([a for b in list(sample_dict.values()) for a in b])]
    df_frac_long[variable] = np.nan
    for i in list(sample_dict.keys()):
        df_frac_long.loc[df_frac_long['Sample'].isin(sample_dict[i]), variable] = i
    df_frac_long['Patient'] = df_frac_long['Sample'].map(sample_pairing_dict)

    if len(cell_type_order) == 0:

        cell_type_order = sorted(df_frac_long['Cell type'].unique())

    if len(variable_order) == 0:

        variable_order = sorted(df_frac_long[variable].unique())

    if not figsize:

        figsize = (len(df_frac_long['Cell type'].unique()), 6)

    if not isinstance(colors, list):

        colors = [colors(i) for i in range(len(variable_order))]
        
    if print_pval:

        text_format='simple'

    else:

        text_format='star'

    df_frac_long.sort_values(['Patient', 'Cell type', variable], inplace=True)

    plt.figure(figsize=figsize)
    ax = sns.stripplot(data=df_frac_long, x='Cell type', y='fraction', hue=variable, jitter=0,
                       order=cell_type_order, hue_order=variable_order, dodge=True, palette=colors)
    
    box_pairs = [((i, variable_order[0]), (i, variable_order[1])) for i in cell_type_order]

    annotator = Annotator(ax=ax, pairs=box_pairs, data=df_frac_long, x='Cell type', y='fraction', 
                          order=cell_type_order, hue=variable, hue_order=variable_order)

    annotator.configure(test=test, comparisons_correction=None, text_format=text_format, loc='inside', show_test_name=False, verbose=0)
    annotator.apply_and_annotate()
    
    plt.xticks(rotation=45, ha='right')
    sns.despine(right=True, top=True)
    plt.xlabel('')
    plt.ylabel('Cell type fraction')

    n = -0.15
    for i in cell_type_order:
        
        x_ = df_frac_long[(df_frac_long['Cell type']==i) & (df_frac_long[variable]==variable_order[0])].sort_values(['Patient', 'Sample'])
        y_ = df_frac_long[(df_frac_long['Cell type']==i) & (df_frac_long[variable]==variable_order[1])].sort_values(['Patient', 'Sample'])
        x = x_[x_['Patient'].isin(list(set(x_['Patient'])&set(y_['Patient'])))]['fraction']
        y = y_[y_['Patient'].isin(list(set(x_['Patient'])&set(y_['Patient'])))]['fraction']

        for j in range(len(x)):
            try:
                plt.plot([n-0.05,n+0.35], [x.tolist()[j], y.tolist()[j]], c='k', linewidth=0.3)
            except:
                continue
        #plt.plot([n-0.15, n+0.1], [y.mean(), y.mean()], c='r', linewidth=2, zorder=10)
        #plt.plot([n+0.2, n+0.45], [x.mean(), x.mean()], c='r', linewidth=2, zorder=10)
        n += 1

def plot_diff_abundance_paired_specific(adata, cell_type_col, cell_type, sample_col, variable, variable_order, 
                                        sample_dict, sample_pairing_dict, test='wilcoxon', 
                                        colors=['lightblue', 'darksalmon'], print_pval=False, one_sided=False, figsize=(4, 5)):
    
    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))
    for i in df[sample_col].unique():
        df_frac[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]
    df_frac = df_frac.divide(df_frac.sum(axis=0), axis=1)

    df_frac_long = df_frac.T
    df_frac_long = df_frac_long.reset_index()
    df_frac_long = pd.melt(df_frac_long, id_vars=['index'], value_vars=df_frac_long.columns[1:])
    df_frac_long = df_frac_long.rename(columns={'index':'Sample', 'variable':'Cell type', 'value':'fraction'})

    df_frac_long = df_frac_long[df_frac_long['Sample'].isin([a for b in list(sample_dict.values()) for a in b])]
    df_frac_long[variable] = np.nan
    for i in list(sample_dict.keys()):
        df_frac_long.loc[df_frac_long['Sample'].isin(sample_dict[i]), variable] = i
    df_frac_long['Patient'] = df_frac_long['Sample'].map(sample_pairing_dict)

    if isinstance(cell_type, str):

        df_frac_long = df_frac_long[df_frac_long['Cell type'] == cell_type]

    elif len(cell_type) == 1:

        df_frac_long_ = df_frac_long.drop_duplicates('Sample').sort_values('Sample')
        
        df_frac_long_['fraction'] = df_frac_long[df_frac_long['Cell type'].isin(cell_type[0])].groupby('Sample').sum()['fraction'].values

        df_frac_long = df_frac_long_

    elif len(cell_type) == 2:

        df_frac_long_ratio = df_frac_long.drop_duplicates('Sample').sort_values('Sample')

        df_frac_long_ratio['fraction'] = df_frac_long[df_frac_long['Cell type'].isin(cell_type[0])].groupby('Sample').sum()['fraction'].values / \
                                          df_frac_long[df_frac_long['Cell type'].isin(cell_type[1])].groupby('Sample').sum()['fraction'].values
        df_frac_long_ratio['fraction'] = np.log2(df_frac_long_ratio['fraction'])

        df_frac_long = df_frac_long_ratio

    plt.figure(figsize=figsize)
    ax = sns.stripplot(data=df_frac_long, x=variable, y='fraction',
                       order=variable_order, dodge=True, jitter=0, palette=colors)
    sns.despine(right=True, top=True)

    x_ = df_frac_long[df_frac_long[variable]==variable_order[0]].sort_values(['Patient', 'Sample'])
    y_ = df_frac_long[df_frac_long[variable]==variable_order[1]].sort_values(['Patient', 'Sample'])
    x = x_[x_['Patient'].isin(list(set(x_['Patient'])&set(y_['Patient'])))]['fraction'].to_numpy()
    y = y_[y_['Patient'].isin(list(set(x_['Patient'])&set(y_['Patient'])))]['fraction'].to_numpy()
    
    x_ = x[(~np.isinf(x))&(~np.isinf(y))]
    y_ = y[(~np.isinf(x))&(~np.isinf(y))]
    x = x_
    y = y_

    if test == 'ttest_rel':
        p = st.ttest_rel(x, y, nan_policy='omit')[1]
    elif test == 'wilcoxon':
        p = st.wilcoxon(x, y, nan_policy='omit')[1]
    if one_sided:
        p = p/2
    if p <= 0.001:
        star = '***'
    elif p <= 0.01:
        star = '**'
    elif p <= 0.05:
        star = '*'
    else:
        star = 'ns'
    if x.mean() < y.mean():
        color='k'
    else:
        color='k'

    for j in range(len(x)):
        plt.plot([0, 1], [x.tolist()[j], y.tolist()[j]], c='k', linewidth=0.3)
        
    x = x[(x<1e4) & (x>-1e4)]
    y = y[(y<1e4) & (y>-1e4)]
    
    val_max = max(max(x), max(y))
    val_min = min(min(x), min(y))
    val_range = val_max - val_min

    if print_pval:
        ax.text(0.5, val_max+val_range*0.02, 'P = %.2g'%p, color=color, ha='center')
    else:
        ax.text(0.5, val_max+val_range*0.02, star, color=color, ha='center')

    plt.xlim(-0.5, 1.5)
    plt.xlabel('')
    plt.ylabel('Cell type fraction')
    
def plot_diff_exp_volcano(df_de_res, level_to_be_contrasted, level_to_contrast, x_metric=None, genes_to_be_marked=None, n_gene=20):
    
    df_de_res = df_de_res.sort_values('pvalue').drop_duplicates('name').set_index('name')
    
    pval_min = df_de_res['-log10pvals'].replace(np.inf, np.nan).dropna().max()
    df_de_res.loc[df_de_res['-log10pvals']==np.inf, '-log10pvals'] = pval_min

    if not x_metric:

        if 'scores' in df_de_res.columns.tolist():
            x_metric = 'scores'
        else:
            x_metric = 'log2FoldChange'

    if x_metric == 'scores':
            x_metric_label = 'DE score'
    elif x_metric == 'log2FoldChange':
            x_metric_label = 'Log2FC'

    log2fc_max = df_de_res[x_metric].abs().max()
    
    if not genes_to_be_marked:
        
        genes_to_be_marked_pos = (st.zscore(df_de_res.query('%s > 0'%x_metric)[x_metric]) * st.zscore(df_de_res.query('%s > 0'%x_metric)['-log10pvals'])).sort_values(ascending=False).index[:n_gene].tolist()
        genes_to_be_marked_neg = (st.zscore(df_de_res.query('%s < 0'%x_metric)[x_metric].abs()) * st.zscore(df_de_res.query('%s < 0'%x_metric)['-log10pvals'])).sort_values(ascending=False).index[:n_gene].tolist()
        genes_to_be_marked = genes_to_be_marked_pos + genes_to_be_marked_neg
    
    plt.figure(figsize=(5, 5))
    
    sns.scatterplot(data=df_de_res, x=x_metric, y='-log10pvals', s=1, linewidth=0, hue='significant', hue_order=['True', 'False'], palette=['darksalmon', 'lightblue'], rasterized=True)
    sns.despine(right=True, top=True)
    
    plt.xlim(-log2fc_max*1.1, log2fc_max*1.1)
    plt.xlabel('%s: %s vs. %s' % (x_metric_label, level_to_be_contrasted, level_to_contrast))
    plt.ylabel('-log10(P-value)')
    
    plt.legend(title='Significant', loc=2)
    
    texts = [plt.text(df_de_res.loc[i, x_metric], \
                      df_de_res.loc[i, '-log10pvals'], \
                      i, fontsize=5) for i in genes_to_be_marked]
    adjust_text(texts, arrowprops=dict(arrowstyle='-', color='k', linewidth=0.5), expand_points=(2, 2))

def plot_pseudobulk_diff_exp_by_cell_type(adata, sample_col, cell_type_col, cell_type_order, gene_order,
                                          test='Wilcoxon', colors=['lightblue', 'darksalmon'], 
                                          figsize=(12, 6), print_pval=False):

    samps = adata.obs[sample_col].unique()
    
    df = pd.DataFrame({})
    
    adata_zscore = st.zscore(adata.raw[:, gene_order].X.toarray(), axis=0, nan_policy='omit')
    
    for gene in gene_order:
    
        df_gene = pd.DataFrame({}, index=samps)

        for cell in cell_type_order:

            df_gene[cell] = [np.mean(adata_zscore[(adata.obs[sample_col]==samp) & (adata.obs[cell_type_col]==cell), gene_order.index(gene)]) for samp in samps]
    
        df_gene['Gene'] = gene
        
        df = pd.concat([df, df_gene], ignore_index=True)
        
    df_long = pd.melt(df.reset_index(), id_vars=['index', 'Gene'], value_vars=df.columns)
    df_long.columns = ['sample', 'gene', 'cell_type', 'expression']
    df_long['expression'] = df_long['expression'].fillna(0)
    
    plt.figure(figsize=figsize)

    ax = sns.stripplot(data=df_long, x='gene', y='expression', hue='cell_type', 
                       order=gene_order, hue_order=cell_type_order, dodge=True, jitter=0, 
                       palette=colors, rasterized=True)

    if print_pval:
    
        box_pairs=[((gene, cell_type_order[0]), (gene, cell_type_order[1]))
                   for gene in df_long['gene'].unique()] + \
                  [((gene, cell_type_order[0]), (gene, cell_type_order[2])) 
                   for gene in df_long['gene'].unique()]

        annotator = Annotator(ax=ax, pairs=box_pairs, data=df_long, x='gene', y='expression', hue='cell_type', 
                            order=gene_order, hue_order=cell_type_order)

        annotator.configure(test=test, comparisons_correction=None, text_format='simple', loc='inside', show_test_name=False, verbose=0)
        annotator.apply_and_annotate()
        
    else:

        box_pairs=[((gene, cell_type_order[0]), (gene, cell_type_order[1]))
                   for gene in df_long['gene'].unique()] + \
                  [((gene, cell_type_order[0]), (gene, cell_type_order[2])) 
                   for gene in df_long['gene'].unique()]

        annotator = Annotator(ax=ax, pairs=box_pairs, data=df_long, x='gene', y='expression', hue='cell_type', 
                            order=gene_order, hue_order=cell_type_order)

        annotator.configure(test=test, comparisons_correction=None, text_format='star', loc='inside', show_test_name=False, verbose=0)
        annotator.apply_and_annotate()
    
    n = 0

    for gene in gene_order:

        if len(cell_type_order) == 3:

            x = df_long.query('gene == @gene & cell_type == @cell_type_order[0]').set_index('sample').loc[samps, 'expression']
            y = df_long.query('gene == @gene & cell_type == @cell_type_order[1]').set_index('sample').loc[samps, 'expression']
            z = df_long.query('gene == @gene & cell_type == @cell_type_order[2]').set_index('sample').loc[samps, 'expression']

            for samp in samps:

                plt.plot([n-0.25, n], [x[samp], y[samp]], c='k', linewidth=0.3)
                plt.plot([n, n+0.25], [y[samp], z[samp]], c='k', linewidth=0.3)

        elif len(cell_type_order) == 2:

            x = df_long.query('gene == @gene & cell_type == @cell_type_order[0]').set_index('sample').loc[samps, 'expression']
            y = df_long.query('gene == @gene & cell_type == @cell_type_order[1]').set_index('sample').loc[samps, 'expression']

            for samp in samps:

                plt.plot([n-0.25, n+0.25], [x.tolist()[samp], y.tolist()[samp]], c='k', linewidth=0.3)

        n += 1

    plt.xticks(rotation=45, ha='right')
    sns.despine(right=True, top=True)
    plt.xlabel('')
    plt.ylabel('Average expression Z-score')

# Visualize VDJ data

def plot_diversity_bar(dict_diversity, metric, sample_subset=[], cell_type_order=None, 
                       variable=None, variable_order=None, variable_dict=None, colors=['salmon'], 
                       test='t-test_welch', box_pairs=None, print_pval=False):

    df_diversity = dict_diversity[metric]
    df_diversity = df_diversity.reset_index().melt(id_vars=['index'], value_vars=df_diversity.columns)
    df_diversity.columns = ['cell_type', 'sample', 'value']

    if len(sample_subset) > 0:
        df_diversity = df_diversity.query('sample in @sample_subset')

    if not variable:

        if not cell_type_order:

            cell_type_order = df_diversity.dropna().groupby('cell_type').count().query('sample >= 2').index
            df_diversity = df_diversity.query('cell_type in @cell_type_order')
            cell_type_order = df_diversity.groupby('cell_type').mean()['value'].sort_values(ascending=False).index

        plt.figure(figsize=(len(cell_type_order), 6))
        sns.barplot(data=df_diversity, x='sample', y='value', 
                    order=cell_type_order, color=colors[0], 
                    errwidth=1, errcolor='k', capsize=0.1)
        plt.xlabel('')
        plt.ylabel(metric.replace('_', ' ').title())
        plt.xticks(rotation=45, ha='right')
        sns.despine(right=True, top=True)

    else:
        
        df_diversity = df_diversity.query('sample in @variable_dict.keys()')

        df_diversity[variable] = df_diversity['sample'].map(variable_dict)

        if not cell_type_order:

            cell_type_order = df_diversity.dropna().groupby('cell_type').count().query('sample >= 2').index
            df_diversity = df_diversity.query('cell_type in @cell_type_order')
            cell_type_order = df_diversity[['cell_type', 'value']].groupby('cell_type').mean()['value'].sort_values(ascending=False).index

        if not variable_order:
            
            variable_order = df_diversity[variable].unique()

        plt.figure(figsize=(len(cell_type_order)*len(set(variable_dict.values()))*0.15, 6))

        ax = sns.barplot(data=df_diversity, x='cell_type', y='value', 
                         hue=variable, order=cell_type_order, hue_order=variable_order, 
                         palette=colors, errwidth=1, errcolor='k', capsize=0.1)

        if box_pairs:

            if print_pval:

                text_format='simple'

            else:

                text_format='star'
                
            box_pairs_ = box_pairs.copy()

            for box_pair in box_pairs_:

                x = df_diversity[df_diversity[variable]==box_pair[0][1]].query('cell_type == @box_pair[0][0]')['value'].dropna()
                y = df_diversity[df_diversity[variable]==box_pair[1][1]].query('cell_type == @box_pair[1][0]')['value'].dropna()

                if not (len(x) >= 2 and len(y) >= 2):
                    
                    box_pairs.remove(box_pair)

            annotator = Annotator(ax=ax, pairs=box_pairs, data=df_diversity, x='cell_type', y='value', 
                                    hue=variable, order=cell_type_order, hue_order=variable_order)

            annotator.configure(test=test, comparisons_correction=None, text_format=text_format, loc='inside', show_test_name=False, verbose=0)
            annotator.apply_and_annotate()

        plt.axhline(y=0, ls='-', lw=2, c='k')
        plt.xlabel('')
        plt.ylabel(metric.replace('_', ' ').title())
        plt.xticks(rotation=45, ha='right')
        plt.legend(title=variable)
        sns.despine(right=True, top=True)

def ci95(x):
    x = x.astype(float)
    x = x[~np.isnan(x)]
    return st.t.interval(alpha=0.95, df=len(x)-1, loc=np.mean(x), scale=st.sem(x))

def plot_differential_expansion_radar(dict_diversity, metric,  
                                      variable, variable_order, variable_dict,
                                      cell_type_order=None, test='Mann-Whitney', colors=['lightblue', 'darksalmon']):

    df_expand = dict_diversity[metric].T.reindex(index=variable_dict[variable_order[0]]+variable_dict[variable_order[1]])

    if not cell_type_order:

        cell_type_order = df_expand.columns[((~df_expand.isna()).sum(axis=0))>=2].tolist()

    df_expand = df_expand[cell_type_order]
    categories = cell_type_order
    
    x = df_expand.reindex(index=variable_dict[variable_order[0]])
    y = df_expand.reindex(index=variable_dict[variable_order[1]])

    N = len(categories)
    angles = [n / float(N) * 2 * np.pi for n in range(N)]
    angles += angles[:1]

    mean = x.mean(axis=0).tolist()
    upper = [ci95(i)[1] for i in x.to_numpy().T]
    lower = [ci95(i)[0] for i in x.to_numpy().T]

    mean += mean[:1]
    upper += upper[:1]
    lower += lower[:1]
    #lower = [0 if i < 0 else i for i in lower]

    mean1 = y.mean(axis=0).tolist()
    upper1 = [ci95(i)[1] for i in y.to_numpy().T]
    lower1 = [ci95(i)[0] for i in y.to_numpy().T]

    mean1 += mean1[:1]
    upper1 += upper1[:1]
    lower1 += lower1[:1]
    #lower1 = [0 if i < 0 else i for i in lower1]
    
    min_val = min(min(np.array(lower)[~np.isnan(lower)&~np.isinf(lower)]), 
                  min(np.array(lower1)[~np.isnan(lower1)&~np.isinf(lower1)]))
    
    max_val = max(max(np.array(upper)[~np.isnan(upper)&~np.isinf(upper)]), 
                  max(np.array(upper1)[~np.isnan(upper1)&~np.isinf(upper1)]))

    plt.figure(figsize=(6, 6))

    plt.subplot(polar=True)
    plt.plot(angles, mean, label='%s-%s'%(variable, variable_order[0]), color=colors[0])
    plt.fill_between(angles, upper, lower, facecolor=colors[0], alpha=0.5)

    ax = plt.subplot(polar=True)
    plt.plot(angles, mean1, label='%s-%s'%(variable, variable_order[1]), color=colors[1])
    plt.fill_between(angles, upper1, lower1, facecolor=colors[1], alpha=0.5)

    ax.set_rlabel_position(0)
    plt.yticks(list(np.arange(min_val, max_val, (max_val-min_val)/4)), 
               ['%.2f'%i for i in list(np.arange(min_val, max_val, (max_val-min_val)/4))], 
               color="grey", size=8, ha='right')

    #ax.spines["start"].set_color("none")
    #ax.spines["polar"].set_color("none")

    ax.yaxis.grid(linestyle='--', linewidth=1)
    ax.xaxis.grid(linestyle='--', linewidth=1)

    plt.xticks(angles[:-1], labels=categories, color='k', size=8, ha='center')

    for cat in categories:

        x_ = x[cat]
        y_ = y[cat]

        if x_.mean() > y_.mean():
            color = colors[0]
        else:
            color = colors[1]

        try:

            if test == 'Mann-Whitney':
                p = st.mannwhitneyu(x_.tolist(), y_.tolist(), nan_policy='omit')[1]
            elif test == 't-test_ind':
                p = st.ttest_ind(x_.tolist(), y_.tolist(), nan_policy='omit')[1]
            elif test == 't-test_welch':
                p = st.ttest_ind(x_.tolist(), y_.tolist(), equal_var=False, nan_policy='omit')[1]
            elif test == 'Wilcoxon':
                p = st.wilcoxon(x_.tolist(), y_.tolist(), nan_policy='omit')[1]
            elif test == 't-test_rel':
                p = st.ttest_rel(x_.tolist(), y_.tolist(), nan_policy='omit')[1]

        except:

            p = 1

        if not p <= 1 or np.isnan(p):

            p = 1

        fontsize=20
        if p > 0.05:
            star = 'ns'
            color = 'k'
            fontsize=10
        elif p < 0.05 and p > 0.01:
            star = '*'
        elif p < 0.01 and p > 0.001:
            star = '**'
        elif p < 0.001:
            star = '***'

        plt.text(angles[categories.index(cat)], 
                 (max_val-min_val)*0.95+min_val,
                 star, color=color, fontsize=fontsize,
                 va='center', ha='center', 
                 rotation=angles[categories.index(cat)]*60-90)

    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))

def plot_repertoire_overlap_heatmap(df_ovlp, method='average', metric='euclidean', order=None, 
                                    vmin=None, vmax=None, cmap=None, linewidth=3, figsize=(9, 9)):   

    if not vmin:

        vmin = df_ovlp.min().min()

    if not vmax:

        vmax = sorted(set(df_ovlp.to_numpy().flatten()))[-2]

    if not cmap:

        cmap = plt.cm.viridis

    if order:

        df_ovlp = df_ovlp.reindex(index=order, columns=order)

        cm = sns.clustermap(df_ovlp.fillna(0), 
                            cmap=cmap, linewidth=linewidth, figsize=figsize, 
                            row_cluster=False, col_cluster=False,
                            vmin=vmin, vmax=vmax
                            )

    else:

        cm = sns.clustermap(df_ovlp.fillna(0), 
                            cmap=cmap, linewidth=linewidth, figsize=figsize, 
                            row_cluster=True, col_cluster=True,
                            vmin=vmin, vmax=vmax
                            )

        cm.ax_row_dendrogram.set_visible(False)
        cm.ax_col_dendrogram.set_visible(False)

    mask = np.tril(np.ones_like(df_ovlp), 0)
    values = cm.ax_heatmap.collections[0].get_array().reshape(df_ovlp.shape)
    new_values = np.ma.array(values, mask=mask)
    new_values = np.ma.array(new_values, mask=df_ovlp.isna())

    cm.ax_heatmap.collections[0].set_array(new_values)
    cm.ax_heatmap.set_xlabel('')
    cm.ax_heatmap.set_ylabel('')
    cm.ax_heatmap.set_xticks([])
    cm.ax_heatmap.set_ylabel('')

    cm.ax_cbar.set_title('Overlap')

def plot_repertoire_overlap_circle(df_ol, vertex_weights=[], vertex_colors=[], figsize=(7, 7), save=None):

    m,n = df_ol.shape

    df_ol[:] = np.where(np.arange(m)[:, None] >= np.arange(n),0,df_ol)

    if len(vertex_weights) == 0:

        vertex_weights = [1]*len(df_ol)

    if len(vertex_colors) == 0:

        vertex_colors = [plt.cm.Set1(i) for i in range(len(df_ol))]

    ro.r("library(CellChat)")

    ro.r("""
            netVisual_circle <-function(net, color.use = NULL,title.name = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, remove.isolate = FALSE, top = 1,
                                        weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = NULL, vertex.label.cex=1,vertex.label.color= "black",
                                        edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                        edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2, vertex.size = NULL,
                                        arrow.width=1,arrow.size = 0.2){
              if (!is.null(vertex.size)) {
                warning("'vertex.size' is deprecated. Use `vertex.weight`")
              }
              if (is.null(vertex.size.max)) {
                if (length(unique(vertex.weight)) == 1) {
                  vertex.size.max <- 5
                } else {
                  vertex.size.max <- 15
                }
              }
              options(warn = -1)
              thresh <- stats::quantile(net, probs = 1-top)
              net[net < thresh] <- 0

              if ((!is.null(sources.use)) | (!is.null(targets.use)) | (!is.null(idents.use)) ) {
                if (is.null(rownames(net))) {
                  stop("The input weighted matrix should have rownames!")
                }
                cells.level <- rownames(net)
                df.net <- reshape2::melt(net, value.name = "value")
                colnames(df.net)[1:2] <- c("source","target")
                # keep the interactions associated with sources and targets of interest
                if (!is.null(sources.use)){
                  if (is.numeric(sources.use)) {
                    sources.use <- cells.level[sources.use]
                  }
                  df.net <- subset(df.net, source %in% sources.use)
                }
                if (!is.null(targets.use)){
                  if (is.numeric(targets.use)) {
                    targets.use <- cells.level[targets.use]
                  }
                  df.net <- subset(df.net, target %in% targets.use)
                }
                if (!is.null(idents.use)) {
                  if (is.numeric(idents.use)) {
                    idents.use <- cells.level[idents.use]
                  }
                  df.net <- filter(df.net, (source %in% idents.use) | (target %in% idents.use))
                }
                df.net$source <- factor(df.net$source, levels = cells.level)
                df.net$target <- factor(df.net$target, levels = cells.level)
                df.net$value[is.na(df.net$value)] <- 0
                net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
              }
              net[is.na(net)] <- 0


              if (remove.isolate) {
                idx1 <- which(Matrix::rowSums(net) == 0)
                idx2 <- which(Matrix::colSums(net) == 0)
                idx <- intersect(idx1, idx2)
                net <- net[-idx, ]
                net <- net[, -idx]
              }

              g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
              edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
              coords<-layout_(g,layout)
              if(nrow(coords)!=1){
                coords_scale=scale(coords)
              }else{
                coords_scale<-coords
              }
              if (is.null(color.use)) {
                color.use = scPalette(length(igraph::V(g)))
              }
              if (is.null(vertex.weight.max)) {
                vertex.weight.max <- max(vertex.weight)
              }
              vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

              loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
              igraph::V(g)$size<-vertex.weight
              igraph::V(g)$color<-color.use[igraph::V(g)]
              igraph::V(g)$frame.color <- color.use[igraph::V(g)]
              igraph::V(g)$label.color <- vertex.label.color
              igraph::V(g)$label.cex<-vertex.label.cex
              if(label.edge){
                igraph::E(g)$label<-igraph::E(g)$weight
                igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
              }
              if (is.null(edge.weight.max)) {
                edge.weight.max <- max(igraph::E(g)$weight)
              }
              if (weight.scale == TRUE) {
                #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
                igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
              }else{
                igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
              }

              igraph::E(g)$arrow.width<-arrow.width
              igraph::E(g)$arrow.size<-arrow.size
              igraph::E(g)$label.color<-edge.label.color
              igraph::E(g)$label.cex<-edge.label.cex
              igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
              igraph::E(g)$loop.angle <- rep(0, length(igraph::E(g)))

              if(sum(edge.start[,2]==edge.start[,1])!=0){
                igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
              }
              radian.rescale <- function(x, start=0, direction=1) {
                c.rotate <- function(x) (x + start) %% (2 * pi) * direction
                c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
              }
              label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
              label.dist <- vertex.weight/max(vertex.weight)+2
              plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
                   vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
              if (!is.null(title.name)) {
                text(0,1.5,title.name, cex = 1.1)
              }
              # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
              # grid.echo()
              # gg <-  grid.grab()
              gg <- recordPlot()
              return(gg)
            }
        """)

    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.globalenv['mat'] = df_ol

    ro.globalenv['vertex_colors'] = pd.Series([matplotlib.colors.to_hex(i) for i in vertex_colors])
    ro.globalenv['vertex_weights'] = pd.Series(vertex_weights)

    if not save:
        
        save = './tmp/trash/tmp.pdf'
        
    ro.globalenv['fig_path'] = save
    ro.globalenv['fig_size'] = pd.Series(figsize)

    ro.r("""
        pdf(fig_path, width=fig_size[0], height=fig_size[1])
        netVisual_circle(as.matrix(mat), color.use = vertex_colors, vertex.weight = vertex_weights, weight.scale = T, 
                         edge.weight.max = max(mat), title.name = NULL, arrow.size=0)
        dev.off()
        """)

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()

def plot_differential_repertoire_overlap_circle(adata, sample_col, cell_type_col, clone_col, 
                                                variable, variable_order, metric='smaller', test='ttest_welch',
                                                sample_subset=[], cell_type_subset=[], sample_dict=None,
                                                figsize=(7, 7), edge_colors=('salmon', 'lightblue'), save=None):

    if sample_dict:

        samps_1 = sample_dict[variable_order[0]]
        samps_2 = sample_dict[variable_order[1]]

    else:

        samps_1 = adata.obs[adata.obs[variable]==variable_order[0]][sample_col].unique().tolist()
        samps_2 = adata.obs[adata.obs[variable]==variable_order[1]][sample_col].unique().tolist()

    df_rep_ol_per_samp_1 = get_repertoire_overlap_per_subset_df(adata, sample_col, cell_type_col, clone_col, 
                                                                sample_subset=samps_1, cell_type_subset=[], metric=metric)
    df_rep_ol_per_samp_2 = get_repertoire_overlap_per_subset_df(adata, sample_col, cell_type_col, clone_col, 
                                                                sample_subset=samps_2, cell_type_subset=[], metric=metric)

    cell_types = list(set(df_rep_ol_per_samp_1['group_1'].tolist()+df_rep_ol_per_samp_1['group_2'].tolist()+df_rep_ol_per_samp_2['group_1'].tolist()+df_rep_ol_per_samp_2['group_2'].tolist()))

    df_rep_ol_diff = pd.DataFrame(0, index=cell_types, columns=cell_types)
    df_rep_ol_pval = pd.DataFrame(1, index=cell_types, columns=cell_types)

    for cell1 in cell_types:
        for cell2 in cell_types:
            try:

                try:
                    x = df_rep_ol_per_samp_1.query('group_1 == @cell1 & group_2 == @cell2').iloc[0, 2:]
                except:
                    x = df_rep_ol_per_samp_1.query('group_1 == @cell2 & group_2 == @cell1').iloc[0, 2:]
                try:
                    y = df_rep_ol_per_samp_2.query('group_1 == @cell1 & group_2 == @cell2').iloc[0, 2:]
                except:
                    y = df_rep_ol_per_samp_2.query('group_1 == @cell2 & group_2 == @cell1').iloc[0, 2:]

                if len(x.dropna()) == 1 or len(y.dropna()) == 1:

                    df_rep_ol_diff.loc[cell1, cell2] = 0
                    df_rep_ol_pval.loc[cell1, cell2] = 1

                else:

                    df_rep_ol_diff.loc[cell1, cell2] = np.mean(x.dropna()) - np.mean(y.dropna())
                    x = x.tolist()
                    y = y.tolist()
                    if test == "ttest_welch":
                        df_rep_ol_pval.loc[cell1, cell2] = st.ttest_ind(x, y, equal_var=False, nan_policy='omit')[1]
                    elif test == "ttest_ind":
                        df_rep_ol_pval.loc[cell1, cell2] = st.ttest_ind(x, y, nan_policy='omit')[1]
                    elif test == "ttest_rel":
                        df_rep_ol_pval.loc[cell1, cell2] = st.ttest_rel(x, y, nan_policy='omit')[1]
                    elif test == "mannwhitney":
                        df_rep_ol_pval.loc[cell1, cell2] = st.mannwhitneyu(x, y, nan_policy='omit')[1]
                    elif test == "wilcoxon":
                        df_rep_ol_pval.loc[cell1, cell2] = st.wilcoxon(x, y, nan_policy='omit')[1]

            except:

                df_rep_ol_diff.loc[cell1, cell2] = 0
                df_rep_ol_pval.loc[cell1, cell2] = 1

    df_rep_ol_res = df_rep_ol_diff*(df_rep_ol_pval<0.05)
    df_rep_ol_res = df_rep_ol_res*np.tril(np.ones_like(df_rep_ol_res), -1)
    df_rep_ol_zero = pd.DataFrame(0, index=df_rep_ol_res.index, columns=df_rep_ol_res.columns)

    ro.r("library(CellChat)")

    ro.r("""
            netVisual_diffInteraction <- function(obj1, obj2, comparison = c(1,2), color.use = NULL, color.edge = c('#b2182b','#2166ac'), title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                                  weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                                                  edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                                  edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
                                                  arrow.width=1,arrow.size = 0.2){
              options(warn = -1)

              net.diff <- obj2 - obj1

              net <- net.diff
              if ((!is.null(sources.use)) | (!is.null(targets.use))) {
                df.net <- reshape2::melt(net, value.name = "value")
                colnames(df.net)[1:2] <- c("source","target")
                # keep the interactions associated with sources and targets of interest
                if (!is.null(sources.use)){
                  if (is.numeric(sources.use)) {
                    sources.use <- rownames(net.diff)[sources.use]
                  }
                  df.net <- subset(df.net, source %in% sources.use)
                }
                if (!is.null(targets.use)){
                  if (is.numeric(targets.use)) {
                    targets.use <- rownames(net.diff)[targets.use]
                  }
                  df.net <- subset(df.net, target %in% targets.use)
                }
                cells.level <- rownames(net.diff)
                df.net$source <- factor(df.net$source, levels = cells.level)
                df.net$target <- factor(df.net$target, levels = cells.level)
                df.net$value[is.na(df.net$value)] <- 0
                net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
                net[is.na(net)] <- 0
              }

              if (remove.isolate) {
                idx1 <- which(Matrix::rowSums(net) == 0)
                idx2 <- which(Matrix::colSums(net) == 0)
                idx <- intersect(idx1, idx2)
                net <- net[-idx, ]
                net <- net[, -idx]
              }

              net[abs(net) < stats::quantile(abs(net), probs = 1-top)] <- 0

              g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
              edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
              coords<-layout_(g,layout)
              if(nrow(coords)!=1){
                coords_scale=scale(coords)
              }else{
                coords_scale<-coords
              }
              if (is.null(color.use)) {
                color.use = scPalette(length(igraph::V(g)))
              }
              if (is.null(vertex.weight.max)) {
                vertex.weight.max <- max(vertex.weight)
              }
              vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5

              loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
              igraph::V(g)$size<-vertex.weight
              igraph::V(g)$color<-color.use[igraph::V(g)]
              igraph::V(g)$frame.color <- color.use[igraph::V(g)]
              igraph::V(g)$label.color <- vertex.label.color
              igraph::V(g)$label.cex<-vertex.label.cex
              if(label.edge){
                igraph::E(g)$label<-igraph::E(g)$weight
                igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
              }
              igraph::E(g)$arrow.width<-arrow.width
              igraph::E(g)$arrow.size<-arrow.size
              igraph::E(g)$label.color<-edge.label.color
              igraph::E(g)$label.cex<-edge.label.cex
              #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
              igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
              igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)

              igraph::E(g)$weight <- abs(igraph::E(g)$weight)

              if (is.null(edge.weight.max)) {
                edge.weight.max <- max(igraph::E(g)$weight)
              }
              if (weight.scale == TRUE) {
                #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
                igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
              }else{
                igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
              }


              if(sum(edge.start[,2]==edge.start[,1])!=0){
                igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
              }
              radian.rescale <- function(x, start=0, direction=1) {
                c.rotate <- function(x) (x + start) %% (2 * pi) * direction
                c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
              }
              label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
              label.dist <- vertex.weight/max(vertex.weight)+2
              plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
                   vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
              if (!is.null(title.name)) {
                text(0,1.5,title.name, cex = 1.1)
              }
              # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
              # grid.echo()
              # gg <-  grid.grab()
              gg <- recordPlot()
              return(gg)
            }
        """)

    anndata2ri.activate()
    anndata2ri.scipy2ri.activate()

    ro.globalenv['mat1'] = df_rep_ol_zero
    ro.globalenv['mat2'] = df_rep_ol_res
    ro.globalenv['vertex_weights'] = adata.obs[adata.obs[sample_col].isin(samps_1+samps_2)].groupby(cell_type_col).count()['total_counts'][cell_types]

    if not save:
        
        save = './tmp/trash/tmp.pdf'
        
    ro.globalenv['fig_path'] = save
    ro.globalenv['fig_size'] = pd.Series(figsize)
    ro.globalenv['edge_colors'] = pd.Series([matplotlib.colors.to_hex(i) for i in edge_colors])
    ro.r("""
        pdf(fig_path, width=fig_size[0], height=fig_size[1])
        netVisual_diffInteraction(as.matrix(mat1), as.matrix(mat2), comparison = c(1,2), vertex.weight = vertex_weights, 
                                weight.scale = T, edge.weight.max = max(abs(as.matrix(mat2)-as.matrix(mat1)))*0.5, 
                                title.name = NULL, arrow.size=0, color.edge = edge_colors)
        dev.off()
        """)

    anndata2ri.scipy2ri.deactivate()
    anndata2ri.deactivate()

def plot_clo_co_exp(adata, sample_col, clone_col, sample,
                    variable, variable_order, num_clo_ol=50, 
                    color='salmon', return_freq=False):

    adata_tis1 = adata.obs[(adata.obs[sample_col]==sample)&(adata.obs[variable]==variable_order[0])].copy()
    adata_tis2 = adata.obs[(adata.obs[sample_col]==sample)&(adata.obs[variable]==variable_order[1])].copy()

    clo_cnt_tis1 = adata_tis1[clone_col].value_counts()
    clo_cnt_tis2 = adata_tis2[clone_col].value_counts()

    clo_cnt_tis1 = clo_cnt_tis1[clo_cnt_tis1>0]
    clo_cnt_tis2 = clo_cnt_tis2[clo_cnt_tis2>0]

    clo_union = set(clo_cnt_tis1.index) | set(clo_cnt_tis2.index)
    clo_ol = list(set(clo_cnt_tis1.index) & set(clo_cnt_tis2.index))

    clo_cnt_tis1 = clo_cnt_tis1.reindex(index=clo_union)
    clo_cnt_tis2 = clo_cnt_tis2.reindex(index=clo_union)

    clo_cnt_tis1 = clo_cnt_tis1[clo_ol]
    clo_cnt_tis2 = clo_cnt_tis2[clo_ol]

    clo_freq_tis1 = clo_cnt_tis1 / clo_cnt_tis1.sum()
    clo_freq_tis2 = clo_cnt_tis2 / clo_cnt_tis2.sum()
    
    if not num_clo_ol:

        num_clo_ol = 0
    
    if len(clo_ol) >= num_clo_ol:

        x = clo_freq_tis1
        y = clo_freq_tis2

        plt.figure(figsize=(5, 5))

        plt.scatter(np.log10(x), np.log10(y), s=(x+y)*2e3, linewidth=1, edgecolor='k', alpha=0.7, color=color)
        sns.regplot(np.log10(x), np.log10(y), scatter_kws={'s':0}, line_kws={'color':'k', 'lw':2})

        vmin = min(min(np.log10(x)), min(np.log10(y)))
        vmax = max(max(np.log10(x)), max(np.log10(y)))
        offset = 0.05*(vmax - vmin)

        plt.plot([vmin-offset, vmax+offset], [vmin-offset, vmax+offset], lw=1, ls='--', color='k')

        plt.xlim(vmin-offset, vmax+offset)
        plt.ylim(vmin-offset, vmax+offset)

        plt.xlabel('%s clone fraction (log 10)'%variable_order[0])
        plt.ylabel('%s clone fraction (log 10)'%variable_order[1])

        plt.text(vmin+offset, vmax-offset, 'Rs = %.2f, p = %.2g'%st.spearmanr(x, y, nan_policy='omit'), size=12)
        
        plt.title(sample)
        
        if return_freq:
            
            return pd.DataFrame({'Freq1':x, 'Freq2':y}, index=clo_ol)
        
    else:
        
        print('%s does not have more than %s ovelapped clones to plot.'%(sample, num_clo_ol))


def plot_composition_complex_log(adata, x, y, x_order=None, y_order=None, x_annot=None, xticklabels=False, \
                             figsize=(10, 8), colors=None, palette=palettable.cartocolors.qualitative.Prism_10.get_mpl_colormap(), \
                             heat_norm=True, bar_x_norm=False, bar_y_norm=True, \
                             bar_width=0.6, bar_x_linkage=True, bar_x_link_alpha=0.5, \
                             heat_size_scale=50, heat_marker='.'):

    # Parse data
    df = adata.obs
    df_frac = pd.DataFrame({}, index=sorted(df[y].unique()))
    if not x_order:
        x_order = list(df[x].unique())
    if not y_order:
        y_order = df_frac.index
    for i in x_order:
        df_frac[i] = df[df[x]==i].groupby(y).count().iloc[:, 0]
    df_frac = df_frac.loc[y_order, :]
    df_frac_heat = df_frac.copy()
    df_frac_bar_x = df_frac.copy()
    df_frac_bar_y = df_frac.copy()
    if heat_norm:
        df_frac_heat = df_frac_heat.divide(df_frac_heat.sum(axis=0), axis=1)*100
        heat_unit = 'proportion'
    else:
        df_frac_heat = df_frac_heat
        heat_unit = 'count'
    if bar_x_norm:
        df_frac_bar_x = np.log10(df_frac_bar_x+1)
        bar_x_ylabel = 'Counts (log10)'
    else:
        df_frac_bar_x = df_frac_bar_x
        bar_x_ylabel = 'Counts'
    if bar_y_norm:
        vec_frac_bar_y = df_frac_bar_y.sum(axis=1)/df_frac_bar_y.sum(axis=0).sum()*100
        bar_y_xlabel = 'Proportion (%)'
    else:
        vec_frac_bar_y = np.log10(df_frac_bar_y.sum(axis=1)+1)
        bar_y_xlabel = 'Count (log10)'
    if not colors:
        palette = [plt.cm.get_cmap(palette)(perc/100) for perc in np.arange(0, 100, 100/len(y_order))]
    else:
        palette = colors
    # Create Fig and gridspec
    fig = plt.figure(figsize=figsize)
    grid = plt.GridSpec(nrows=4, ncols=4, hspace=0.05, wspace=0.05)
    # Define the axes
    ax_heatmap_main = fig.add_subplot(grid[1:3, 1:3])
    ax_heatmap_main_legend = fig.add_subplot(grid[3, 0])
    if x_annot: 0
    ax_heatmap_x_annot = fig.add_subplot(grid[3, 1:3])
    ax_heatmap_y_annot = fig.add_subplot(grid[1:3, 0])
    ax_bar_x = fig.add_subplot(grid[0, 1:3])
#    ax_bar_y = fig.add_subplot(grid[1:3, 3])
    # Configure the axes
    ax_heatmap_main.set_zorder(10)
    ax_heatmap_main.tick_params(axis='y', which='major', pad=15)
    #ax_heatmap_main.add_patch(Rectangle((-0.5, -0.5), 8, 8,
    #                  alpha=1, fill=None, linewidth=2))
    #ax_heatmap_y_annot.margins(y=0.01)
    ax_heatmap_y_annot.set_ylim(-0.5, len(df_frac_bar_y.index)-0.5)
    ax_heatmap_y_annot.axis('off')
    #ax_bar_x.margins(x=0.01)
    ax_bar_x.set_ylabel(bar_x_ylabel)
    ax_bar_x.set_xlim(-0.5, len(df_frac_bar_x.columns)-0.5)
    ax_bar_x.get_xaxis().set_visible(False)
    ax_bar_x.locator_params(nbins=5, axis='y')
    sns.despine(ax=ax_bar_x, left=False, right=True, top=True, bottom=True)
    ###############
    # Draw the main heatmap
    fancy_heatmap(fig=fig, ax=(ax_heatmap_main, ax_heatmap_main_legend), df=df_frac_heat, only_size=True, unit=heat_unit, \
                  x_order=x_order, y_order=y_order, \
                  size_scale=heat_size_scale, palette=palette, marker=heat_marker)
    if not xticklabels:
        ax_heatmap_main.set_xticks([])
    ################
    # Draw the bottom annotation heatmap
    if x_annot:
        df_x_annot = adata.obs.drop_duplicates(x).set_index(x).loc[x_order, x_annot]
        for i in x_annot:
            if isinstance(df_x_annot[i][0], str):
                df_x_annot[i] = df_x_annot[i].map(dict(zip(df_x_annot[i].unique(), range(len(df_x_annot[i].unique()))))).tolist()
        df_x_annot = df_x_annot.T
        df_x_annot = st.zscore(df_x_annot, axis=1)
        sns.heatmap(df_x_annot, ax=ax_heatmap_x_annot, linewidth=0.5, cbar=False, cmap=plt.cm.Spectral_r)
        ax_heatmap_main.get_xaxis().set_visible(False)
        ax_heatmap_x_annot.set_xlabel('')
        ax_heatmap_x_annot.yaxis.set_label_position("right")
        ax_heatmap_x_annot.yaxis.tick_right()
        if xticklabels:
            ax_heatmap_x_annot.set_xticks(list(range(len(x_order))), x_order, rotation=45, ha='right', va='top')
        else:
            ax_heatmap_x_annot.set_xticks([])
        ax_heatmap_x_annot.set_yticks([i+0.5 for i in list(range(len(x_annot)))], x_annot, rotation=0, ha='left', va='center')
    ################
    # Draw the top stacked bar plot
    fancy_barplot(ax=ax_bar_x, df=df_frac_bar_x, legend=False, \
                  bar_cmap=palette, bar_width=bar_width, bar_linkage=bar_x_linkage, bar_link_alpha=bar_x_link_alpha)

##############################