import sys
import os
sys.path.append(os.path.abspath('./OV_utility_functions'))
from basic_imports import *

import scirpy as ir

from itertools import product, combinations

from scipy.stats.mstats import gmean
from scipy.spatial.distance import cosine, squareform

import igraph as ig
from scipy import sparse
from scipy.sparse import csr_matrix, spmatrix

from tcrdist.repertoire import TCRrep
from clustcr import Clustering

# Bulk

def combine_vdj_data(vdj_dir, cols_to_use, cols_renamed):

	globbed_files = glob.glob(vdj_dir+'/*')

	data = []

	for f in globbed_files:
		dataframe = pd.read_table(f)
		dataframe['sample'] = os.path.basename(f).split('.')[0]
		try: 
			int(dataframe['sample'][0])
			dataframe['sample'] = 'S'+dataframe['sample']
		except:
			dataframe['sample'] = dataframe['sample']
		data.append(dataframe)

	combined_data = pd.concat(data)[['sample']+cols_to_use]
	combined_data.columns = ['sample']+cols_renamed

	return combined_data

def compute_basics(df):

	df_ = df.copy()

	df_count = df_.groupby(['sample']).agg(
		{'#count': 'sum'}).reset_index().rename(columns={'#count': "reads_count"})

	df_diversity = df_.groupby(['sample'],
							  sort=False).size().reset_index(name='clonotype_count')

	df_mean_frequency = df_.groupby(['sample']).agg(
		{'freq': 'mean'}).reset_index().rename(columns={'freq': "mean_frequency"})

	samples = df_['sample'].unique()
	df_geomean_frequency = pd.DataFrame(columns=['sample', 'geomean_frequency'])
	for sample in samples:
		tmp = df_[df_['sample'] == sample]
		geomean_frequency = gmean(tmp['freq'])
		df_data = pd.DataFrame({'sample': sample, 'geomean_frequency': geomean_frequency}, index=[0])
		df_geomean_frequency = pd.concat([df_geomean_frequency, df_data], copy=False, ignore_index=True)

	df_['length_weighted'] = df_['cdr3nt'].str.len()*df_['freq']
	df_mean_cdr3nt_length = df_.groupby(['sample']).agg(
		{'length_weighted': 'sum'}).reset_index().rename(columns={'length_weighted': "mean_cdr3nt_length"})

	# Count unique CDR3
	df_unique_CDR3 = df_.groupby(['cdr3aa', 'sample'], as_index=False)[
		'cdr3nt'].agg({'count': 'count'})

	# Calculate the mean of the unique CDR3 count in each sample
	df_unique_CDR3_mean = df_unique_CDR3.groupby(['sample']).agg(
		{'count': 'mean'}).reset_index().rename(columns={'count': "convergence"})

	# CDR3 nucleotide length
	df_['nt_length'] = df_['cdr3nt'].str.len()

	# Calculate spectratype
	df_spectratype = df_.groupby(['sample', 'nt_length']).agg(
		{'freq': 'sum'}).reset_index().rename(columns={'freq': "spectratype"})

	df_spectratype

	# Merge df_count and df_geomean_frequency first
	df_geomean_frequency = df_geomean_frequency.merge(
		df_count, on='sample', how='left')

	# Create a dataframe that combines all the basic analysis (except for spectratype)
	dfs = [df_diversity, df_mean_frequency, df_geomean_frequency,
		   df_mean_cdr3nt_length, df_unique_CDR3_mean]

	axis = ['sample']

	df_basics = pd.merge(dfs[0], dfs[1], left_on=axis, right_on=axis, how='outer')
	for d in dfs[2:]:
		df_basics = pd.merge(df_basics, d, left_on=axis, right_on=axis, how='outer')

	return df_basics

def cluster_tcr_leiden(df, dist_cutoff= 50, chunk_size = 100, res = 1, n_iter = -1):

	def standardize_gene_name(i):

		if '-' in i:
			if len(i.split('-')[-1]) == 1:
				i = i.split('-')[0] + '-0' + i.split('-')[-1]

		if len(i.split('-')[0].split('R')[1]) == 3:
			try:
				i = i.split('R')[0]+'R'+i.split('R')[1].split('-')[0][:2]+'0'+i.split('R')[1].split('-')[0][-1]+'-'+i.split('-')[1]
			except:
				i = i.split('R')[0]+'R'+i.split('R')[1].split('-')[0][:2]+'0'+i.split('R')[1].split('-')[0][-1]

		return i

	def distance_to_connectivity(distances):

		connectivities = distances.copy()
		d = connectivities.data
		d[d < 0] = 0
		max_value = np.max(d) + 1

		# structure of the matrix stays the same, we can safely change the data only
		connectivities.data = (max_value - d) / max_value
		connectivities.eliminate_zeros()

		return connectivities

	def get_igraph_from_adjacency(adj):

		sources, targets = adj.nonzero()
		weights = adj[sources, targets]

		if isinstance(weights, np.matrix):
			weights = weights.A1
		if isinstance(weights, csr_matrix):
			# this is the case when len(sources) == len(targets) == 0, see #236
			weights = weights.toarray()

		g = ig.Graph(directed=False)
		g.add_vertices(adj.shape[0])  # this adds adjacency.shape[0] vertices
		g.add_edges(list(zip(sources, targets)))

		g.es["weight"] = weights

		return g

	df_ = df.copy()
	df_.columns = ['subject', 'freq', 'count', 'cdr3_b_aa', 'cdr3_b_nucseq', 'v_b_gene', 'd_b_gene', 'j_b_gene']

	df_['v_b_gene'] = [i.replace('TR', 'TCR') if 'TCR' not in i else i for i in df_['v_b_gene']]
	df_['j_b_gene'] = [i.replace('TR', 'TCR') if 'TCR' not in i else i for i in df_['j_b_gene']]

	df_['v_b_gene'] = df_['v_b_gene'].apply(standardize_gene_name)
	df_['j_b_gene'] = df_['j_b_gene'].apply(standardize_gene_name)

	conv = pd.read_csv('~/anaconda3/lib/python3.8/site-packages/tcrdist/db/adaptive_imgt_mapping.csv')
	conv_dict = conv.set_index('adaptive')['imgt'].to_dict()

	df_['v_b_gene'] = df_['v_b_gene'].map(conv_dict)
	df_['j_b_gene'] = df_['j_b_gene'].map(conv_dict)

	df_.dropna(inplace=True)

	tr = TCRrep(cell_df = df_, 
				organism = 'mouse', 
				chains = ['beta'], 
				db_file = 'alphabeta_gammadelta_db.tsv',
				cpus=1,
				compute_distances = False
				)

	tr.cpus = 32
	tr.compute_sparse_rect_distances(radius = dist_cutoff, chunk_size = chunk_size)

	graph = get_igraph_from_adjacency(distance_to_connectivity(tr.rw_beta))

	leiden = graph.community_leiden(objective_function="modularity",
									 resolution_parameter=res,
									 n_iterations=n_iter)

	tr.cell_df['id'] = tr.cell_df['cdr3_b_aa'] + tr.cell_df['v_b_gene'] + tr.cell_df['j_b_gene']

	tr.clone_df['id'] = tr.clone_df['cdr3_b_aa'] + tr.clone_df['v_b_gene'] + tr.clone_df['j_b_gene']
	tr.clone_df['leiden_cluster'] = leiden.membership

	tr.cell_df['leiden_cluster'] = tr.cell_df['id'].map(tr.clone_df.set_index('id')['leiden_cluster'].to_dict())

	tr.cell_df['new_id'] = tr.cell_df['subject'] + 'C' + tr.cell_df['leiden_cluster'].astype(str)
	tr.cell_df['freq'] = tr.cell_df['new_id'].map(tr.cell_df[['new_id', 'freq']].groupby('new_id').sum()['freq'].to_dict())
	tr.cell_df['count'] = tr.cell_df['new_id'].map(tr.cell_df[['new_id', 'count']].groupby('new_id').sum()['count'].to_dict())
	tr.cell_df.drop_duplicates('new_id', inplace=True)

	df_comb_clustered = tr.cell_df[['subject', 'freq', 'count', 'leiden_cluster']]
	df_comb_clustered.columns = ['sample', 'freq', '#count', 'cluster']
	df_comb_clustered.dropna(inplace=True)
	df_comb_clustered['freq'] = (df_comb_clustered.set_index('sample')['freq']/df_comb_clustered.groupby('sample').sum().loc[df_comb_clustered['sample'], 'freq']).tolist()
	df_comb_clustered['cluster'] = df_comb_clustered['cluster'].astype(int).astype(str)

	return df_comb_clustered

def cluster_tcr_faiss_mcl(df, method='mcl', n_cpus=8):
	
	clustering = Clustering(method=method, n_cpus=n_cpus)

	cdr3 = df['cdr3aa']
	cdr3.columns = ['junction_aa']
	output = clustering.fit(cdr3)

	output.clusters_df
	
	df_comb_clustered = df.copy()
	df_comb_clustered['clustcr_cluster'] = df_comb_clustered['cdr3aa'].map(output.clusters_df.set_index('junction_aa')['cluster'].to_dict())
	
	df_comb_clustered['new_id'] = df_comb_clustered['sample'] + 'C' + df_comb_clustered['clustcr_cluster'].astype(str)
	df_comb_clustered['freq'] = df_comb_clustered['new_id'].map(df_comb_clustered[['new_id', 'freq']].groupby('new_id').sum()['freq'].to_dict())
	df_comb_clustered['#count'] = df_comb_clustered['new_id'].map(df_comb_clustered[['new_id', '#count']].groupby('new_id').sum()['#count'].to_dict())
	df_comb_clustered.drop_duplicates('new_id', inplace=True)

	df_comb_clustered = df_comb_clustered[['sample', 'freq', '#count', 'clustcr_cluster']]
	df_comb_clustered = df_comb_clustered.rename(columns={'clustcr_cluster':'cluster'})
	df_comb_clustered.dropna(inplace=True)
	df_comb_clustered['freq'] = (df_comb_clustered.set_index('sample')['freq']/df_comb_clustered.groupby('sample').sum().loc[df_comb_clustered['sample'], 'freq']).tolist()
	df_comb_clustered['cluster'] = df_comb_clustered['cluster'].astype(int).astype(str)
	return df_comb_clustered

def compute_diversity_metrics(df):

	df_diversity = df.groupby(['sample'], sort=False).size().reset_index(name='clonotype_count')
	df_ = pd.merge(df, df_diversity, on=['sample'])

	# Shannon-Wiener index
	# Calculation step 1
	df_['shannon_index'] = -(df_['freq']*np.log(df_['freq']))

	# Calculation step 2
	df_shannon = df_.groupby(['sample']).agg(
		{'shannon_index': 'sum'}).reset_index().rename(columns={'': "shannon_index"})

	# Calculation step 3, Shannon-Wienex index is shown in the shannon_wiener_index column
	df_shannon['shannon_wiener_index'] = np.exp(df_shannon['shannon_index'])
	df_shannon_index = df_shannon[['sample', 'shannon_wiener_index']]

	# Normalized Shannon-Wiener index
	# Calculation step 1 - merge df_shannon and df_diversity (which contains clonotype counts)
	df_shannon = pd.merge(df_shannon, df_diversity, on=[
						  'sample'])

	# Calculation step 2 - calculate normalized Shannon-Wienex index, it is shown in the normalized_shannon_wiener_index column
	df_shannon['normalized_shannon_wiener_index'] = df_shannon['shannon_index'] / \
		np.log(df_shannon['clonotype_count'])

	df_norm_shannon = df_shannon[['sample','shannon_wiener_index', 'normalized_shannon_wiener_index']]

	# Inverse Simpson index
	# Calculation step 1
	df_['simpson_index'] = (df_['freq']**2)

	# Calculation step 2
	df_simpson = df_.groupby(['sample']).agg(
		{'simpson_index': 'sum'}).reset_index().rename(columns={'': "simpson_index"})

	# Calculation step 3, Inverse Simpson index is shown in the inverse_simpson_index column
	df_simpson['inverse_simpson_index'] = 1/df_simpson['simpson_index']

	# Gini Simpson index
	df_simpson['gini_simpson_index'] = 1-df_simpson['simpson_index']

	# D50 index
	# Create an empty dataframe for storing results
	df_D50 = pd.DataFrame()

	# Create a list of the sample names
	samples = set(df_['sample'])

	for sample in samples:

		# Store the rows related to the sample
		df_temp = df_.loc[df_['sample'] == sample]

		# Sort the sample clonotypes by frequency in descending order
		df_temp = df_temp.sort_values(by='freq', ascending=False)

		# Create a column to store the order
		df_temp['clonotype_number'] = np.arange(df_temp.shape[0])+1

		# Compute and store the cumulative sum of the frequencies
		df_temp['accum_freq'] = df_temp['freq'].cumsum()

		# Find out the first accumulated frequency that is above 50%
		df_temp = df_temp.loc[(df_temp['accum_freq'] >= 0.5)
							  & (df_temp['accum_freq'] <= 0.6)]
		df_temp = df_temp.head(1)

		# Calculate D50 index and store in the result dataframe
		df_temp = df_temp.head(1)
		df_temp['D50_index'] = df_temp['clonotype_number'] / \
			df_temp['clonotype_count']*100
		df_D50 = pd.concat([df_D50, df_temp])

	df_D50 = df_D50[['sample', 'D50_index']]

	# chao1 and chao1_SD
	# Create an empty dataframe for storing results
	df_chao1 = pd.DataFrame()

	# Get the columns needed for calculation from df
	df1 = df_[['sample', '#count', 'clonotype_count']]

	# Create a list of the sample names
	samples = set(df1['sample'])

	for sample in samples:

		# Store the rows related to the sample
		df_temp = df1.loc[df1['sample'] == sample]

		try:

			# Count singleton in the sample
			singleton = len(df_temp.loc[df_temp['#count'] == 1])

			# Count doubleton in the sample
			doubleton = len(df_temp.loc[df_temp['#count'] == 2])

			# Calculate Chao1 estimate
			chao1 = int(df_temp['clonotype_count'].values[0]) + \
				((singleton * (singleton-1))/(2*(doubleton+1)))
			df_temp['chao1'] = chao1

			# Calculate Chao1 estimate standard deviation
			step1 = 1/4*((singleton/doubleton)**4)
			step2 = (singleton/doubleton)**3
			step3 = 1/2*((singleton/doubleton)**2)
			step4 = doubleton * (step1+step2+step3)
			df_temp['chao1_SD'] = step4**(1/2)

		except:

			df_temp['chao1'] = np.nan
			df_temp['chao1_SD'] = np.nan

		# Store the results in the result dataframe
		df_chao1 = pd.concat([df_chao1, df_temp], axis=0, sort=False)

		# Remove the duplicates results in the result dataframe
		df_chao1 = df_chao1[['sample', 'chao1', 'chao1_SD']]
		df_chao1 = df_chao1.drop_duplicates(subset=['sample'], keep='first')

	# Gini coefficient
	# Create an empty dataframe for storing results
	df_gini = pd.DataFrame()

	# Create a list of the sample names
	samples = set(df['sample'])

	for sample in samples:

		# Store the rows related to the sample
		df_temp = df.loc[df['sample'] == sample]

		def gini(list_of_values):
			sorted_list = sorted(list_of_values)
			height, area = 0, 0
			for value in sorted_list:
				height += value
				area += height - value / 2.
			fair_area = height * len(list_of_values) / 2.
			return (fair_area - area) / fair_area

		# Calculate gini coefficient
		df_temp['gini_coefficient'] = gini(df_temp['freq'])

		# Store the results in the result dataframe
		df_gini = pd.concat([df_gini, df_temp], sort=False)

		# Remove the duplicates results in the result dataframe
		df_gini = df_gini[['sample', 'gini_coefficient']]
		df_gini = df_gini.drop_duplicates(subset=['sample'], keep='first')

	# Clonality
	# Add the clonotype counts as a column to the dataframe
	df_clonocount = df.groupby(['sample'], sort=False).size().reset_index(name='clonotype_count')
	df_clonality = pd.merge(df, df_clonocount, on=['sample'])

	# Calculate 1-Pielou index
	df_clonality['clonality'] = df_clonality['freq']*np.log(df_clonality['freq'])/np.log(df_clonality['clonotype_count'])
	df_clonality= df_clonality.groupby(['sample']).agg({'clonality':'sum'}).reset_index().rename(columns={'':"clonality"})
	df_clonality['1_pielou'] = df_clonality['clonality'] + 1
	df_clonality = df_clonality[['sample', 'clonality','1_pielou']]

	# Combine all metrics
	# Create a dataframe that combines all the diversity analysis
	dfs = [df_norm_shannon, df_simpson, df_D50, df_chao1, df_gini, df_clonality]

	df_diversity = pd.merge(dfs[0], dfs[1], left_on=['sample'], right_on=[
						   'sample'], how='outer')

	for d in dfs[2:]:
		df_diversity = pd.merge(df_diversity, d, left_on=['sample'], right_on=[
							   'sample'], how='outer')

	df_diversity.set_index('sample', inplace=True)

	return df_diversity

def compute_relative_abundance(df, plot=False, sample_order=[], cutoffs=[0.01, 0.001, 0.0001, 0.00001]):

	df_ = df.copy()

	# Define clonotype groups based on frequency
	def clonotype_group (row):
		if row['freq'] > cutoffs[0] and row['freq'] <= 1:
			return 'Hyperexpanded'
		if row['freq'] > cutoffs[1] and row['freq'] <= cutoffs[0]: 
			return 'Large'
		if row['freq'] > cutoffs[2] and row['freq'] <= cutoffs[1]:
			return 'Medium'
		if row['freq'] > cutoffs[3] and row['freq'] <= cutoffs[2]:
			return 'Small'
		if row['freq'] > 0 and row['freq'] <= cutoffs[3]:
			return 'Rare'

	# Apply the clonotype_group function to the dataframe
	df_['clonotype_group'] = df_.apply(lambda row: clonotype_group(row), axis=1)

	# Calculate the relative abundance in each sample based on clonotype groups
	df_relative_abundance= df_.groupby(['sample','clonotype_group']).agg({'freq':'sum'}).reset_index().rename(columns={'':"relative_abundance"})

	if plot:

		label_order = ['Hyperexpanded', 'Large', 'Medium', 'Small', 'Rare']
		if len(sample_order) == 0:
			sample_order = sorted(df_relative_abundance['sample'].unique())
		
		colors = plt.cm.Oranges_r(np.linspace(0, 1, 5))
		ax = df_relative_abundance.groupby(['sample','clonotype_group'])['freq'].sum().unstack().loc[sample_order, label_order].plot(kind='bar', stacked=True, color=colors)
		ax.set_xlabel('sample',fontsize=20)
		ax.set_ylabel('clonotype frequency',fontsize=20)
		plt.xticks(fontsize=20)
		plt.yticks(fontsize=20)
		plt.gcf().set_size_inches(20,10)
		sns.despine()

	else:

		return df_relative_abundance

def compute_repertoire_overlap(df):

	df_ = df.copy()

	if 'cluster' in df_.columns:
		df_compare = pd.merge(df_, df_, on=['cluster'], suffixes=['_1', '_2'])
	else:
		# If there are more than one V or J gene, leave only the first one
		df_['v'].str.replace("(,).*", "", regex=True)
		df_['j'].str.replace("(,).*", "", regex=True)
		df_compare = pd.merge(df_, df_, on=['cdr3nt', 'cdr3aa', 'v', 'd', 'j'], suffixes=['_1', '_2'])

	df_compare = df_compare[df_compare['sample_1'] != df_compare['sample_2']]

	df_compare['#count_1**2'] = df_compare['#count_1'] * df_compare['#count_1']
	df_compare['#count_2**2'] = df_compare['#count_2'] * df_compare['#count_2']
	df_compare['#count_1*2'] = df_compare['#count_1'] * df_compare['#count_2']

	# Prepare data for future analysis
	sample_names = df_["sample"].drop_duplicates()
	sample_names_size = sample_names.size
	sample_names.index = np.arange(0, sample_names_size)

	samples = {}

	for sample in sample_names:
		filtered_samples = df_.loc[(df_['sample'] == sample)]
		samples[sample] = filtered_samples

	df_overlaps = {}

	for i in range(0, sample_names_size):
		sample1 = sample_names[i]

		df_compare1 = (df_compare['sample_1'] == sample1)

		for j in range(i + 1, sample_names_size):
			sample2 = sample_names[j]

			df_compare2 = (df_compare['sample_2'] == sample2)

			df_rows = df_compare.loc[df_compare1 & df_compare2]
			df_overlaps[f'{sample1}:{sample2}'] = df_rows

	# Create an empty dataframe for storing results
	columns = ['sample_1', 'sample_2', 
			   'jaccard_index', 'overlap_coefficient', 'morisita_horn_index', 'tversky_index', 'cosine_similarity', 
			   'pearson_correlation_count', 'pearson_correlation_freq', 'relative_overlap_diversity', 
			   'geometric_mean_of_relative_overlap_frequencies', 'clonotype_wise_sum_of_geometric_mean_frequencies', 
			   'jensen_shannon_divergence_v_usage']

	df_ovlp_metrics = pd.DataFrame(columns=columns)

	for i in range(0, sample_names_size):

		sample1 = sample_names[i]
		df_sample1 = samples[sample1]

		data = []

		for j in range(i + 1, sample_names_size):

			sample2 = sample_names[j]
			df_sample2 = samples[sample2]

			# Calculate Jaccard Index
			jaccard_overlap = df_overlaps[f'{sample1}:{sample2}'].shape[0]
			jaccard = jaccard_overlap / \
				(df_sample1.shape[0] + df_sample2.shape[0] - jaccard_overlap)

			# Calculate Overlap Coefficient
			sample_overlap = df_overlaps[f'{sample1}:{sample2}']
			overlap_coefficient = sample_overlap.shape[0] / \
				min(df_sample1.shape[0], df_sample2.shape[0])

			# Calculate Morisita-Horn index
			sum_sample_1_count = df_sample1["#count"].sum()
			sum_sample_2_count = df_sample2["#count"].sum()

			sample_overlap = df_overlaps[f'{sample1}:{sample2}']

			sum_sample_1 = sample_overlap["#count_1**2"].sum()
			sum_sample_2 = sample_overlap["#count_2**2"].sum()

			sum_count_product = sample_overlap["#count_1*2"].sum()
			step1 = (sum_sample_1 / ((sum_sample_1_count)**2)) + \
				(sum_sample_2 / ((sum_sample_2_count)**2))
			step2 = step1 * sum_sample_1_count * sum_sample_2_count
			step3 = 2 * sum_count_product
			morisita_horn_index = step3 / step2

			# Calculate Tversky index
			df_overlap = df_overlaps[f'{sample1}:{sample2}']

			tversky_index = df_overlap.shape[0] / (df_overlap.shape[0] + 0.5 * (
				df_sample1.shape[0] - df_overlap.shape[0]) + 0.5 * (df_sample2.shape[0] - df_overlap.shape[0]))		 

			# Calculate cosine similarity
			cos_vec = (1 - cosine(df_overlap["freq_1"], df_overlap["freq_2"])) 

			# Calculate Pearson correlation on clonotype count
			pearson_correlation_cnt = df_overlap['#count_1'].corr(df_overlap['#count_2'])

			# Calculate Pearson correlation on clonotype frequency
			pearson_correlation_freq = df_overlap['freq_1'].corr(df_overlap['freq_2'])

			# Calculate relative overlap diversity
			relative_overlap_diversity = df_overlap.shape[0] / (df_sample1.shape[0] * df_sample2.shape[0])

			# Calculate geo mean of rel ovlp freq
			geo_mean_rel_ovlp_freq = (df_overlap['freq_1'].sum() * df_overlap['freq_2'].sum())**0.5

			# Calculate clono-wise sum of geo mean of freq
			clono_sum_geo_mean_freq = ((df_overlap['freq_1']*df_overlap['freq_2'])**0.5).sum()

			# Calculate Jensen-Shannon divergence of v usage
			try:

				def kl_divergence(p, q):
					return -np.sum(p * np.log2(q / p))

				def js_divergence(p, q):
					m = (1 / 2) * (p + q)
					return (1 / 2) * kl_divergence(p, m) + (1 / 2) * kl_divergence(q, m)

				df_sample1_v = df_sample1.groupby(['v'], as_index=False)[
					'freq'].agg({'sumfreq_1': 'sum'})
				df_sample2_v = df_sample2.groupby(['v'], as_index=False)[
					'freq'].agg({'sumfreq_2': 'sum'})

				df_JSD_combine = pd.merge(df_sample1_v, df_sample2_v, on=['v'])

				JSD_sample1 = df_JSD_combine[['sumfreq_1']].to_numpy()
				JSD_sample2 = df_JSD_combine[['sumfreq_2']].to_numpy()

				js_div_v_use = js_divergence(JSD_sample1, JSD_sample2)

			except:

				js_div_v_use = np.nan

			data.append({'sample_1': sample1, 'sample_2': sample2,
						'jaccard_index':jaccard, 
						'overlap_coefficient':overlap_coefficient, 
						'morisita_horn_index':morisita_horn_index, 
						'tversky_index':tversky_index, 
						'cosine_similarity':cos_vec, 
						'pearson_correlation_count':pearson_correlation_cnt, 
						'pearson_correlation_freq':pearson_correlation_freq, 
						'relative_overlap_diversity':relative_overlap_diversity, 
						'geometric_mean_of_relative_overlap_frequencies':geo_mean_rel_ovlp_freq, 
						'clonotype_wise_sum_of_geometric_mean_frequencies':clono_sum_geo_mean_freq, 
						'jensen_shannon_divergence_v_usage':js_div_v_use})

			data.append({'sample_1': sample2, 'sample_2': sample1,
						'jaccard_index':jaccard, 
						'overlap_coefficient':overlap_coefficient, 
						'morisita_horn_index':morisita_horn_index, 
						'tversky_index':tversky_index, 
						'cosine_similarity':cos_vec, 
						'pearson_correlation_count':pearson_correlation_cnt, 
						'pearson_correlation_freq':pearson_correlation_freq, 
						'relative_overlap_diversity':relative_overlap_diversity, 
						'geometric_mean_of_relative_overlap_frequencies':geo_mean_rel_ovlp_freq, 
						'clonotype_wise_sum_of_geometric_mean_frequencies':clono_sum_geo_mean_freq, 
						'jensen_shannon_divergence_v_usage':js_div_v_use})

		df_data = pd.DataFrame(data)
		df_ovlp_metrics = pd.concat([df_ovlp_metrics, df_data], copy=False, ignore_index=True)

	return df_ovlp_metrics

def compute_motifs(df, aa_len=8, nt_len=24):

	df_ = df.copy()

	# Add CDR3 amino acid length as a new column
	df_['aa_length'] = df_['cdr3aa'].str.len()

	df_aa_spectratype = df_.groupby(['sample', 'aa_length'], as_index=False)['freq'].agg({'spectratype': 'sum'})
	df_aa_max_spectratype = df_aa_spectratype.loc[df_aa_spectratype.groupby('sample')['spectratype'].idxmax()]

	# Define the function to count amino acid motifs (k is the length of the motif)
	def aamotif(k, aa_list):
		aamotifCount = {}
		for aa in aa_list:
			for i in range(len(aa)-k+1):
				aamotif = aa[i:i+k]
				aamotifCount[aamotif] = aamotifCount.get(aamotif, 0)+1
		return aamotifCount

	# Create an empty dataframe for storing results
	df_aa_motif = pd.DataFrame()

	# Create a list of the sample names
	samples = set(df_['sample'])

	for sample in samples:

		# Store the rows related to the sample
		df_temp = df_.loc[df_['sample'] == sample]

		# Use amino acid motif length of 6 as an example
		df_temp = aamotif(aa_len, df_temp['cdr3aa'])
		df_temp = pd.DataFrame(df_temp.items(), columns=['motif', 'count'])
		df_temp['sample'] = sample

		# Append the dataframe based on amino acid motifs and stores in the result dataframe
		df_aa_motif = df_aa_motif.append(df_temp, ignore_index=True)

		# Add the hospitalization information as a column
		df_aa_motif_1 = df_aa_motif.merge(df_aa_max_spectratype[['sample']], on='sample')

	# Define the function to count nucleotide motifs (k is the length of the motif)
	def ntmotif(k, nt_list):
		ntmotifCount = {}
		for nt in nt_list:
			for i in range(len(nt)-k+1):
				ntmotif = nt[i:i+k]
				ntmotifCount[ntmotif] = ntmotifCount.get(ntmotif, 0)+1
		return ntmotifCount

	# Create an empty dataframe for storing results
	df_nt_motif = pd.DataFrame()

	# Create a list of the sample names
	samples = set(df_['sample'])

	for sample in samples:

		# Store the rows related to the sample
		df_temp = df_.loc[df_['sample'] == sample]

		# Use amino acid motif length of 6 as an example
		df_temp = ntmotif(nt_len, df_temp['cdr3nt'])
		df_temp = pd.DataFrame(df_temp.items(), columns=['motif', 'count'])
		df_temp['sample'] = sample

		# Append the dataframe based on amino acid motifs and stores in the result dataframe
		df_nt_motif = df_nt_motif.append(df_temp, ignore_index=True)

		# Add the hospitalization information as a column
		df_nt_motif_1 = df_nt_motif.merge(
			df_aa_max_spectratype[['sample']], on='sample')

	return df_nt_motif_1, df_aa_motif_1

# Single-cell

def gini_coef(list_of_values):
	sorted_list = sorted(list_of_values)
	height, area = 0, 0
	for value in sorted_list:
		height += value
		area += height - value / 2.
	fair_area = height * len(list_of_values) / 2.
	return (fair_area - area) / fair_area

def get_diversity_df(adata, sample_col, cell_type_col, clone_col, filtering=False, cell_num_thresh=50):

	adata_ = adata.copy()

	if filtering:

		df = adata_.obs
		df = df[~df[clone_col].isna()]

		df_cnt = pd.DataFrame({}, index=sorted(df[cell_type_col].unique()))

		for i in df[sample_col].unique():

			df_cnt[i] = df[df[sample_col]==i].groupby(cell_type_col).count().iloc[:, 0]

	dict_diversity = {}

	samples = adata_.obs[sample_col].unique()

	for metric in ['chao1', 'D50', 'gini_index', 'pielou_e', 'simpson_e', 'gini_coef'] + \
				  ['doublet_multiplet', 'multiplet']:
		
		dict_diversity[metric] = pd.DataFrame({}, 
											  index=adata_.obs[cell_type_col].unique(),
											  columns=samples)
		
		if not metric.endswith('let') and not metric == "gini_coef":

			for sample in samples:
				
				try:

					dict_diversity[metric][sample] = ir.tl.alpha_diversity(adata_[adata_.obs[sample_col] == sample], 
																		   groupby=cell_type_col, target_col=clone_col, 
																		   metric=metric, inplace=False)
				except:
					
					continue

		elif metric == 'gini_coef':

			for sample in samples:
				
				try:

					dict_diversity[metric][sample] = ir.tl.alpha_diversity(adata_[adata_.obs[sample_col] == sample], 
																		   groupby=cell_type_col, target_col=clone_col, 
																		   metric=gini_coef, inplace=False)
				except:
					
					continue
				
		elif metric == 'doublet_multiplet':
			
			for sample in samples:
				
				try:
				
					dict_diversity[metric][sample] = ir.tl.summarize_clonal_expansion(adata_[adata_.obs[sample_col] == sample], cell_type_col, 
																					  target_col=clone_col, 
																					  summarize_by='clone_id', normalize=True).iloc[:, 1:].sum(axis=1)
				except:
					
					continue
					
		elif metric == 'multiplet':
			
			for sample in samples:
				
				try:
				
					dict_diversity[metric][sample] = ir.tl.summarize_clonal_expansion(adata_[adata_.obs[sample_col] == sample], cell_type_col, 
																					  target_col=clone_col, 
																					  summarize_by='clone_id', normalize=True).iloc[:, 2:].sum(axis=1)
				except:
					
					continue

		if filtering:

			dict_diversity[metric] = dict_diversity[metric].mask(df_cnt <= cell_num_thresh, np.nan)

		dict_diversity[metric] = dict_diversity[metric].astype(float)

	return dict_diversity

def calc_ovlp(x, y, metric='jaccard'):

	x = [i for i in x if isinstance(i, str) or isinstance(i, int)]
	y = [i for i in y if isinstance(i, str) or isinstance(i, int)]

	if metric == 'jaccard':

		try:

			ovlp = len(set(x)&set(y)) / len(set(x)|set(y))

		except:

			ovlp = np.nan

	elif metric == 'smaller':

		try:

			ovlp = len(set(x)&set(y)) / min(len(set(x)), len(set(y)))

		except:

			ovlp = np.nan

	return ovlp

def get_repertoire_overlap_df(adata, sample_col, cell_type_col, clone_col, 
							  sample_subset=[], cell_type_subset=[], metric='smaller'):

	df_clone = adata.obs[[sample_col, cell_type_col, clone_col]].astype(str)

	if len(sample_subset) > 0:

		df_clone = df_clone[df_clone[sample_col].isin(sample_subset)]

	if len(cell_type_subset) > 0:

		df_clone = df_clone[df_clone[cell_type_col].isin(cell_type_subset)]

	#df_clone = df_clone.loc[:, ~df_clone.columns.duplicated()]

	groups = df_clone[cell_type_col].unique()

	df_ovlp = pd.DataFrame({}, index=groups, columns=groups)

	for group1 in groups:
		
		for group2 in groups:

			x = df_clone[df_clone[cell_type_col]==group1][clone_col]
			y = df_clone[df_clone[cell_type_col]==group2][clone_col]
			df_ovlp.loc[group1, group2] = calc_ovlp(x, y, metric)

	df_ovlp = df_ovlp.astype(float)

	return df_ovlp

def get_repertoire_overlap_per_subset_df(adata, sample_col, cell_type_col, clone_col, 
										 sample_subset=[], cell_type_subset=[], metric='smaller'):

	adata_ = adata.copy()

	if len(sample_subset) == 0:

		sample_subset = adata_.obs[sample_col].unique().tolist()

	else:

		adata_ = adata_[adata_.obs[sample_col].isin(sample_subset)]

	if len(cell_type_subset) == 0:

		cell_type_subset = adata_.obs[cell_type_col].unique().tolist()

	else:

		adata_ = adata_[adata_.obs[cell_type_col].isin(cell_type_subset)]
	
	group_comb = list(combinations(cell_type_subset, 2))

	df_rep_ol_per_subset = pd.DataFrame({'group_1':[i[0] for i in group_comb], 
								  		 'group_2':[i[1] for i in group_comb]},
								  		 columns=['group_1', 'group_2'] + sample_subset) 

	for subset in sample_subset:

		try:

			df_rep_ol_subset = get_repertoire_overlap_df(adata_, sample_col, cell_type_col, clone_col,
														 sample_subset=[subset], cell_type_subset=[], metric=metric)

			df_rep_ol_per_subset[subset] = [df_rep_ol_subset.reindex(index=[i[0]], columns=[i[1]]).iloc[0, 0] for i in group_comb]

		except:

			df_rep_ol_per_subset[subset] = np.nan

	df_rep_ol_per_subset = df_rep_ol_per_subset.iloc[:, :2].join(df_rep_ol_per_subset.iloc[:, 2:].astype(float))
			
	return df_rep_ol_per_subset

def get_repertoire_overlap_combo_df(adata, sample_col, clone_col, 
                                    variable_levels, variable_values, 
                                    sample_subset=[], metric='jaccard'):

    adata_ = adata.copy()
    if 'group' in adata_.obs.columns.tolist():
        del adata_.obs['group']

    for variable_level in variable_levels:
        adata_ = adata_[adata_.obs[variable_level].isin([i[variable_levels.index(variable_level)] for i in variable_values])]
        if 'group' not in adata_.obs.columns.tolist():
            adata_.obs['group'] = adata_.obs[variable_level].astype(str)
        else:
            adata_.obs['group'] = adata_.obs['group'].astype(str) + ' - ' + adata_.obs[variable_level].astype(str)

    df_ol_all = get_repertoire_overlap_per_subset_df(adata_, sample_col, 'group', clone_col, 
                                                     sample_subset=sample_subset, cell_type_subset=[' - '.join(i) for i in variable_values], metric=metric)

    df_ol_all.index = [df_ol_all.loc[i, 'group_1']+' : '+df_ol_all.loc[i, 'group_2'] 
                       for i in df_ol_all.index]

    df_ol_all = df_ol_all.loc[df_ol_all.iloc[:, 2:].mean(axis=1).sort_values(ascending=False).index, ]
    df_ol_all_long = df_ol_all.reset_index().melt(id_vars=['index'], value_vars=df_ol_all.columns[2:])
    df_ol_all_long.columns = ['overlap', 'sample', 'value']
    df_ol_all_long.overlap = \
    pd.Categorical(df_ol_all_long.overlap, 
                   categories=df_ol_all_long.groupby('overlap').mean()['value'].sort_values(ascending=False).index.tolist())
    df_ol_all_long = df_ol_all_long.sort_values('overlap')

    df_ol_mean = df_ol_all.copy()
    df_ol_mean['mean'] = df_ol_mean.iloc[:, 2:].mean(axis=1)
    df_ol_mean = df_ol_mean[['group_1', 'group_2', 'mean']]

    df_ol_mean = \
    pd.concat([df_ol_mean, 
               pd.DataFrame([[i, i, 1] 
                             for i in 
                             set(df_ol_mean['group_1'].tolist()+df_ol_mean['group_2'].tolist())],
                             columns=['group_1', 'group_2', 'mean'])
              ])

    df_ol_mean = df_ol_mean.pivot(index='group_1', columns='group_2', values='mean')

    for i in df_ol_mean.index:
        for j in df_ol_mean.columns:
            if not df_ol_mean.loc[i, j] <= 1:
                df_ol_mean.loc[i, j] = df_ol_mean.loc[j, i]

    df_ol_mean.fillna(0, inplace=True)
    
    return df_ol_all, df_ol_all_long, df_ol_mean

# Calculate Treg score
def calc_tirp_score(adata, trb_v_gene_col, trb_cdr3_aa_col):

	df_trb = adata.obs[[trb_v_gene_col, trb_cdr3_aa_col]]
	df_trb['cell'] = df_trb.index

	anndata2ri.activate()

	ro.r("""
	source('/priv18data1/hliang1_group/yluo5/bins/TiRP/results/TiRP.R')
	""")

	ro.globalenv['input'] = df_trb

	ro.r("""
	res = TiRP(input)
	""")

	adata.obs['TiRP'] = ro.r['res'][['cell', 'TiRP']].set_index('cell').reindex(index=df_trb.index)['TiRP'].tolist()

	anndata2ri.deactivate()