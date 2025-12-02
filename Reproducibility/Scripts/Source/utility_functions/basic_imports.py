import sys
import os
import glob
import pickle
import ntpath
import subprocess
import warnings
warnings.filterwarnings("ignore")
import logging
logging.getLogger('fontTools.subset').setLevel(logging.WARNING)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
from matplotlib import rcParams

import seaborn as sns
import palettable
from adjustText import adjust_text

import pandas as pd
import scipy.stats as st
import numpy as np

from matplotlib.offsetbox import AnchoredText
def add_at(ax, label, size=6, loc=2):
    fp = dict(size=size)
    _at = AnchoredText(label, loc=loc, prop=fp)
    ax.add_artist(_at)
    return _at

from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool

# gencode43 = pd.read_table('./genome/hg38/gencode/gencode.v43.annotation.gtf', skiprows=5, header=None)
# gencode43 = gencode43[gencode43[2]=='transcript']
# gencode43['tx'] = [i.split('transcript_id "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43['gene'] = [i.split('gene_id "')[-1].split('"')[0].split('.')[0] for i in gencode43[8]]
# gencode43['name'] = [i.split('gene_name "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43['gene_type'] = [i.split('gene_type "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43['tx_type'] = [i.split('transcript_type "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43['status'] = [i.split('transcript_status "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43['support'] = [i.split('transcript_support_level "')[-1].split('"')[0] for i in gencode43[8]]
# gencode43.drop_duplicates('tx', inplace=True)
# gencode43.head()

# symbol_conv = pd.read_table('./Miscellaneous/gene_list/human_gene_all_alias_previous_hgnc.txt').fillna('None')
# symbol_conv['Approved symbol'] = symbol_conv['Approved symbol'].apply(lambda x: x if 'orf' in x else x.upper())
# symbol_conv['Alias symbol'] = symbol_conv['Alias symbol'].apply(lambda x: x if 'orf' in x else x.upper())
# symbol_conv['Previous symbol'] = symbol_conv['Previous symbol'].apply(lambda x: x if 'orf' in x else x.upper())
# symbol_conv.head()

def conv_symbol(symbols):
    symbols_keep = set(symbols) & set(symbol_conv['Approved symbol'])
    symbols_to_change = set(symbols) - set(symbol_conv['Approved symbol'])
    symbols_previous_to_change = set(symbols_to_change) & set(symbol_conv['Previous symbol']) 
    symbols_alias_to_change = set(symbols_to_change) & set(symbol_conv['Alias symbol']) 
    df_previous = symbol_conv.drop_duplicates('Previous symbol').set_index('Previous symbol').loc[symbols_previous_to_change, 'Approved symbol']
    df_previous = df_previous.loc[[i not in symbols_keep for i in df_previous], ]
    df_alias = symbol_conv.drop_duplicates('Alias symbol').set_index('Alias symbol').loc[symbols_alias_to_change, 'Approved symbol']
    df_alias = df_alias.loc[[i not in symbols_keep for i in df_alias], ]
    dict_previous = df_previous.to_dict()
    dict_alias = df_alias.to_dict()
    symbols_changed = [dict_previous[i] if i in list(dict_previous.keys()) else dict_alias[i] if i in list(dict_alias.keys()) else i for i in symbols]
    return symbols_changed

def multipletests_omit_nan(arr):
    df_arr = pd.DataFrame(arr)
    non_nan_indices = df_arr.dropna().index
    fdr = multipletests(df_arr.dropna().iloc[:, 0], alpha=0.05, method='fdr_bh')[1]
    arr_new = []
    n = 0
    for i in range(len(arr)):
        if i not in non_nan_indices:
            arr_new.append(np.nan)
        else:
            arr_new.append(fdr[n])
            n += 1
    return arr_new

def corr_omit_nan(x, y):
    x = np.array(x)
    y = np.array(y)
    nas = np.logical_or(np.isnan(x), np.isnan(y))
    try:
        x = x[~nas]
        y = y[~nas]
        corr = st.pearsonr(x, y)
        return corr
    except:
        return [np.nan, np.nan]

def omit_nan(x, y):
    x = np.array(x)
    y = np.array(y)
    nas = np.logical_or(np.isnan(x), np.isnan(y))
    x = x[~nas]
    y = y[~nas]
    return x, y
