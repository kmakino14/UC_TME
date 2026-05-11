setwd(path_to_wd)
source("Reproducibility/Scripts/Source/my.source.R")
source("Reproducibility/Scripts/Source/Seurat_source.R")
source("Reproducibility/Scripts/Source/MultiNicheNet_plot_source.R")
options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(viridis))
suppressMessages(library(cowplot))
suppressMessages(library(dplyr))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggthemes))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(nichenetr))
suppressMessages(library(multinichenetr))
set.seed(1234)

################################################################################################################################
# Step 0: Preparation of the analysis
################################################################################################################################
data_dir = "Reproducibility/Data"

lr_network_all = 
    readRDS("~/reference/multinichenetr/lr_network_human_allInfo_30112033.rds") %>% 
    mutate(
      ligand = convert_alias_to_symbols(ligand, organism = "human"), 
      receptor = convert_alias_to_symbols(receptor, organism = "human"))

lr_network_all = lr_network_all  %>% 
    mutate(ligand = make.names(ligand), receptor = make.names(receptor)) 

lr_network = lr_network_all %>% 
    distinct(ligand, receptor)

ligand_target_matrix = readRDS("~/reference/multinichenetr/ligand_target_matrix_nsga2r_final.rds")

lineage = 'BCG'

DOGMA = file.path(data_dir, "Seurat", paste0("UC_DOGMA_seurat_obj_", lineage, ".rds")) %>% readRDS()

DefaultAssay(DOGMA) = 'RNA'
scRNA = DietSeurat(DOGMA, assay = 'RNA')

Idents(scRNA) = "sample"
BCG = subset(scRNA, ident = c("BC_011","BC_039","BC_023","BC_044","BC_033","BC_048","BC_037","BC_047")) 

BCG$timepoint = fct_collapse(BCG$sample, 
                             pre  = c("BC_011","BC_023","BC_033","BC_037"),
                             post = c("BC_039","BC_044","BC_048","BC_047")
                             )

BCG$patient = fct_collapse(BCG$sample, 
                           P1 = c("BC_011","BC_039"),
                           P2 = c("BC_023","BC_044"),
                           P3 = c("BC_033","BC_048"),
                           P4 = c("BC_037","BC_047")
                           )

BCG$sample_id = paste0(BCG$patient,BCG$timepoint)

sce = as.SingleCellExperiment(BCG)
sce = alias_to_symbol_SCE(sce, "human") %>% makenames_SCE()

################################################################################################################################
# Step 1: Prepare the cell-cell communication analysis
################################################################################################################################

sample_id = "sample_id"       
group_id = "timepoint"
celltype_id = "coarse_celltype"
covariates = "patient"
batches = NA

SummarizedExperiment::colData(sce)$coarse_celltype = SummarizedExperiment::colData(sce)$coarse_celltype %>% make.names()
senders_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()
receivers_oi = SummarizedExperiment::colData(sce)[,celltype_id] %>% unique()

# Define the sender and receiver cell types of interest.
contrasts_oi = c("'post-pre','pre-post'")
contrast_tbl = tibble(contrast = c("post-pre", "pre-post"),
                      group = c("post", "pre"))
sce = sce[, SummarizedExperiment::colData(sce)[,group_id] %in% contrast_tbl$group]

#################################################################################
## Cell-type filtering: determine which cell types are sufficiently present
#################################################################################

min_cells = 10  # each cluster min cell number basically 10, if not permitted set to 5.

abundance_info = get_abundance_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  batches = batches
  )

abundance_info$abund_plot_sample

# Cell type filtering based on cell type abundance information
sample_group_celltype_df = abundance_info$abundance_data %>% 
  dplyr::filter(n > min_cells) %>% 
  ungroup() %>% 
  distinct(sample_id, group_id) %>% 
  cross_join(
    abundance_info$abundance_data %>% 
      ungroup() %>% 
      distinct(celltype_id)
    ) %>% 
  arrange(sample_id)

abundance_df = sample_group_celltype_df %>% left_join(
  abundance_info$abundance_data %>% ungroup()
  )

abundance_df$n[is.na(abundance_df$n)] = 0
abundance_df$keep[is.na(abundance_df$keep)] = FALSE
abundance_df_summarized = abundance_df %>% 
  mutate(keep = as.logical(keep)) %>% 
  group_by(group_id, celltype_id) %>% 
  summarise(samples_present = sum((keep)))

celltypes_absent_one_condition = abundance_df_summarized %>% 
  dplyr::filter(samples_present == 0) %>% pull(celltype_id) %>% unique() 
# find truly condition-specific cell types by searching for cell types 
# truely absent in at least one condition

celltypes_present_one_condition = abundance_df_summarized %>% 
  dplyr::filter(samples_present >= 2) %>% pull(celltype_id) %>% unique() 
# require presence in at least 2 samples of one group so 
# it is really present in at least one condition

condition_specific_celltypes = intersect(
  celltypes_absent_one_condition, 
  celltypes_present_one_condition)

total_nr_conditions = SummarizedExperiment::colData(sce)[,group_id] %>% 
  unique() %>% length() 

absent_celltypes = abundance_df_summarized %>% 
  dplyr::filter(samples_present < 2) %>% 
  group_by(celltype_id) %>% 
  dplyr::count() %>% 
  dplyr::filter(n == total_nr_conditions) %>% 
  pull(celltype_id)
  
print("condition-specific celltypes:")
print(condition_specific_celltypes)
  
print("absent celltypes:")
print(absent_celltypes)

analyse_condition_specific_celltypes = FALSE

if(analyse_condition_specific_celltypes == TRUE){
  senders_oi = senders_oi %>% setdiff(absent_celltypes)
  receivers_oi = receivers_oi %>% setdiff(absent_celltypes)
} else {
  senders_oi = senders_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
  receivers_oi = receivers_oi %>% 
    setdiff(union(absent_celltypes, condition_specific_celltypes))
}
sce = sce[, SummarizedExperiment::colData(sce)[,celltype_id] %in% 
            c(senders_oi, receivers_oi)
          ]

#################################################################################################
## Gene filtering: determine which genes are sufficiently expressed in each present cell type
#################################################################################################

min_sample_prop = 0.50
fraction_cutoff = 0.05

frq_list = get_frac_exprs(
  sce = sce, 
  sample_id = sample_id, celltype_id =  celltype_id, group_id = group_id, 
  batches = batches, 
  min_cells = min_cells, 
  fraction_cutoff = fraction_cutoff, min_sample_prop = min_sample_prop)

genes_oi = frq_list$expressed_df %>% 
  dplyr::filter(expressed == TRUE) %>% pull(gene) %>% unique() 
sce = sce[genes_oi, ]

#################################################################################################
## Pseudobulk expression calculation
#################################################################################################

abundance_expression_info = process_abundance_expression_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  min_cells = min_cells, 
  senders_oi = senders_oi, receivers_oi = receivers_oi, 
  lr_network = lr_network, 
  batches = batches, 
  frq_list = frq_list, 
  abundance_info = abundance_info)

abundance_expression_info$celltype_info$pb_df %>% head()
abundance_expression_info$celltype_info$pb_df_group %>% head()
abundance_expression_info$sender_receiver_info$pb_df %>% head()
abundance_expression_info$sender_receiver_info$pb_df_group %>% head()

#################################################################################################
## Differential expression (DE) analysis
#################################################################################################

DE_info = get_DE_info(
  sce = sce, 
  sample_id = sample_id, group_id = group_id, celltype_id = celltype_id, 
  batches = batches, covariates = covariates, 
  contrasts_oi = contrasts_oi, 
  min_cells = min_cells, 
  expressed_df = frq_list$expressed_df)

# Check DE results
DE_info$celltype_de$de_output_tidy %>% head()
DE_info$hist_pvals

# If you do not want to change p-val cutoff, set to empirical_pval==FALSE
empirical_pval = FALSE
if(empirical_pval == TRUE){
  DE_info_emp = get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
  celltype_de = DE_info_emp$de_output_tidy_emp %>% select(-p_val, -p_adj) %>% 
    rename(p_val = p_emp, p_adj = p_adj_emp)
} else {
  celltype_de = DE_info$celltype_de$de_output_tidy
} 

#################################################################################################
## Combine DE information for ligand-senders and receptors-receivers
#################################################################################################
sender_receiver_de = combine_sender_receiver_de(
  sender_de = celltype_de,
  receiver_de = celltype_de,
  senders_oi = senders_oi,
  receivers_oi = receivers_oi,
  lr_network = lr_network
)
sender_receiver_de %>% head(20)

#################################################################################################
## Ligand activity prediction
#################################################################################################
logFC_threshold = 0.50

p_val_threshold = 0.05
p_val_adj = FALSE 

geneset_assessment = contrast_tbl$contrast %>% 
  lapply(
    process_geneset_data, 
    celltype_de, logFC_threshold, p_val_adj, p_val_threshold
  ) %>% 
  bind_rows() 

print(geneset_assessment, width = Inf)

#################################################################################################
## Perform the ligand activity analysis and ligand-target inference
#################################################################################################
top_n_target = 250
verbose = TRUE
cores_system = 8
n.cores = min(cores_system, celltype_de$cluster_id %>% unique() %>% length()) 

ligand_activities_targets_DEgenes = suppressMessages(suppressWarnings(
  get_ligand_activities_targets_DEgenes(
    receiver_de = celltype_de,
    receivers_oi = intersect(receivers_oi, celltype_de$cluster_id %>% unique()),
    ligand_target_matrix = ligand_target_matrix,
    logFC_threshold = logFC_threshold,
    p_val_threshold = p_val_threshold,
    p_val_adj = p_val_adj,
    top_n_target = top_n_target,
    verbose = verbose, 
    n.cores = n.cores
  )
))

ligand_activities_targets_DEgenes$ligand_activities %>% head(20)

#################################################################################################
## Prioritization: rank cell-cell communication patterns through multi-criteria prioritization
#################################################################################################
ligand_activity_down = FALSE
sender_receiver_tbl = sender_receiver_de %>% distinct(sender, receiver)

metadata_combined = SummarizedExperiment::colData(sce) %>% tibble::as_tibble()

if(!is.na(batches)){
  grouping_tbl = metadata_combined[,c(sample_id, group_id, batches)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group",batches)
} else {
  grouping_tbl = metadata_combined[,c(sample_id, group_id)] %>% 
    tibble::as_tibble() %>% distinct()
  colnames(grouping_tbl) = c("sample","group")
}

prioritization_tables = suppressMessages(generate_prioritization_tables(
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de = sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    contrast_tbl = contrast_tbl,
    sender_receiver_tbl = sender_receiver_tbl,
    grouping_tbl = grouping_tbl,
    scenario = "regular", # all prioritization criteria will be weighted equally
    fraction_cutoff = fraction_cutoff, 
    abundance_data_receiver = abundance_expression_info$abundance_data_receiver,
    abundance_data_sender = abundance_expression_info$abundance_data_sender,
    ligand_activity_down = ligand_activity_down
  ))

# Check the output tables
prioritization_tables$group_prioritization_tbl %>% head(20)

#################################################################################################
## Calculate the across-samples expression correlation between ligand-receptor pairs and target genes
#################################################################################################
lr_target_prior_cor = lr_target_prior_cor_inference(
  receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(), 
  abundance_expression_info = abundance_expression_info, 
  celltype_de = celltype_de, 
  grouping_tbl = grouping_tbl, 
  prioritization_tables = prioritization_tables, 
  ligand_target_matrix = ligand_target_matrix, 
  logFC_threshold = logFC_threshold, 
  p_val_threshold = p_val_threshold, 
  p_val_adj = p_val_adj
  )

#################################################################################################
## Save all the output of MultiNicheNet
#################################################################################################
path = "Reproducibility/Results/MultiNicheNet/"

multinichenet_output = list(
    celltype_info = abundance_expression_info$celltype_info,
    celltype_de = celltype_de,
    sender_receiver_info = abundance_expression_info$sender_receiver_info,
    sender_receiver_de =  sender_receiver_de,
    ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
    prioritization_tables = prioritization_tables,
    grouping_tbl = grouping_tbl,
    lr_target_prior_cor = lr_target_prior_cor
  ) 
multinichenet_output = make_lite_output(multinichenet_output)

save = TRUE
if(save == TRUE){
  saveRDS(multinichenet_output, paste0(path, "multinichenet_output.rds"))
}

## Load
multinichenet_output = paste0(path, "multinichenet_output.rds") %>% readRDS()

contrasts_oi = c("'post-pre','pre-post'")
contrast_tbl = tibble(contrast =
                        c("post-pre", "pre-post"),
                      group = c("post", "pre"))

celltype_order = c("CD4_Tconv","Treg","CD8_T","NK_ILC","B","Mono_Mac","DC")
color_order = c("#0D31D1","#C18500","#799D00","#00AB6E","#00A9BE","#D55E00","#E16A86")
names(color_order) = celltype_order

#==================================================#
## Interpreting the MultiNicheNet analysis output ##
#==================================================#
#################################################################################################
## Interpretable bubble plots
#################################################################################################
## set receivers_oi
group_oi = "post"

prioritized_tbl_oi_Tumor_100 = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  100, 
  groups_oi = group_oi, 
  receivers_oi = c("CD4_Tconv",'Mono_Mac'))

# Select important interactions
prioritized_tbl_oi_Tumor_100 = 
  dplyr::filter(prioritized_tbl_oi_Tumor_100, !ligand %in% c('TGFB1','IL7','CLEC2B','CLEC2B','CSF1','CD48','CD86','TNF','COL5A1','A2M','CEACAM1',
                'ITGAL','FURIN','CLEC2D','CD28','PTPRC','IL16','HMGB1','MMP9','C1QB','AREG','CD40','IL18','CXCL10','F13A1','RARRES2','SERPING1',
                'TNFSF14','C3','APOE','ICAM3','CD59','PZP','CCL4L2')) %>%
  dplyr::filter(., !id %in% c("CP_SLC40A1_Treg_Mono_Mac"))

plot_oi = make_sample_lr_prod_activity_plots_Omnipath2(
  multinichenet_output$prioritization_tables, 
  prioritized_tbl_oi_Tumor_100 %>% inner_join(lr_network_all))

paste0("Reproducibility/Results/Plots/BCG/FigureS12E.pdf") %>% pdf(., w=15,h=18)
 plot_oi
dev.off()

#################################################################################################
## Visualization of differential ligand-target links
#################################################################################################
source("~/MultiNicheNet_plot_source.R")

group_oi = "post"
receiver_oi = "CD4_Tconv"

# specify interested ligands
ligands_oi = c("IL15",'IL23A',"CXCL9",'CXCL10','CXCL16','CCL20','CCL5','CD86','CD28','CD244',
               'ICAM1','VCAM1','BTLA','CD72','ITGB2')

prioritized_tbl_ligands_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  50,
  groups_oi = group_oi, 
  receivers_oi = receiver_oi
  ) %>% dplyr::filter(ligand %in% ligands_oi) # ligands should still be in the output tables of course

# Ligand-Target heatmap
ligand_activities_targets_DEgenes = multinichenet_output$ligand_activities_targets_DEgenes
best_upstream_ligands = prioritized_tbl_ligands_oi$ligand %>% unique()
active_ligand_target_links_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% 
        dplyr::inner_join(contrast_tbl) %>% 
        dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% 
        dplyr::ungroup() %>% dplyr::select(ligand, target, ligand_target_weight, direction_regulation) %>% 
        dplyr::rename(weight = ligand_target_weight )

active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
  cutoff = 0
} else {
  cutoff = 0.2
}

active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, 
                                                                            ligand_target_matrix = ligand_target_matrix, 
                                                                            cutoff = cutoff)

#order_targets_ = active_ligand_target_links_df$target %>% unique() %>% 
#                 generics::intersect(rownames(active_ligand_target_links))
#order_targets = order_targets_ %>% make.names()

order_targets <- c(
  "IFNG","TNF","TNFSF14","IL18R1","IL12RB1","IL23R","IL4I1",
  "GZMA","GZMH",
  "CD40LG","SLAMF1","CD2","CD38","DPP4",
  "CCL5","CCL20","CXCR3","CXCR6","CCR5","CCR6","ITGAL",
  "TBX21","STAT1"
)

ligand_target_matrix2 = ligand_target_matrix
ligand_target_matrix2[ligand_target_matrix2 > 0.5*max(ligand_target_matrix2)] <- 0.5*max(ligand_target_matrix2)

combined_plot = make_ligand_activity_target_plot2(
  group_oi, 
  rev(ligands_oi),
  order_targets,
  receiver_oi, 
  prioritized_tbl_ligands_oi, 
  multinichenet_output$prioritization_tables, 
  multinichenet_output$ligand_activities_targets_DEgenes, 
  contrast_tbl, 
  multinichenet_output$grouping_tbl, 
  multinichenet_output$celltype_info, 
  ligand_target_matrix2, 
  plot_legend = FALSE)

combined_plot$combined_plot

paste0("Reproducibility/Results/Plots/BCG/Figure6G.pdf") %>% pdf(., w=15,h=8)
 combined_plot
dev.off()

#################################################################################################
## Intercellular regulatory network inference and visualization
#################################################################################################

genes = c('CD44','SPP1','VCAN','THBS1','IL1A','IL1B','IFITM1','CD40LG','CXCL10','CD86','IFNG','IL16','TNF','CD244','CXCR6','CXCL16',
  'IL18','IL15','CCL20','ITGB2','CXCL9','CCL5','CD244','CXCR3','TNFRSF14','IL23A','IL18RAP','IL12RB2','IL23R','IL18R1','DPP4','CD28','CCL22'
  )

prioritized_tbl_oi = get_top_n_lr_pairs(
  multinichenet_output$prioritization_tables, 
  100, 
  rank_per_group = FALSE) %>% dplyr::filter(ligand %in% genes | receptor %in% genes) 

lr_target_prior_cor_filtered = 
  multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>% unique() %>% 
  lapply(function(group_oi){
    lr_target_prior_cor_filtered = multinichenet_output$lr_target_prior_cor %>%
      inner_join(
        multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
          distinct(ligand, target, direction_regulation, contrast)
        ) %>% 
      inner_join(contrast_tbl) %>% dplyr::filter(group == group_oi)
    
    lr_target_prior_cor_filtered_up = lr_target_prior_cor_filtered %>% 
      dplyr::filter(direction_regulation == "up") %>% 
      dplyr::filter( (rank_of_target < top_n_target) & (pearson > 0.30 | spearman > 0.30))
    
    lr_target_prior_cor_filtered_down = lr_target_prior_cor_filtered %>% 
      dplyr::filter(direction_regulation == "down") %>% 
      dplyr::filter( (rank_of_target < top_n_target) & (pearson < -0.30 | spearman < -0.30))

    lr_target_prior_cor_filtered = bind_rows(
      lr_target_prior_cor_filtered_up, 
      lr_target_prior_cor_filtered_down
      )
}) %>% bind_rows()

lr_target_df = lr_target_prior_cor_filtered %>% 
  distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation)
network = infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi)
network$links %>% head()
network$nodes %>% head()

network_graph = visualize_network(network, color_order)

paste0("Reproducibility/Results/Plots/BCG/Figure6H.pdf") %>% pdf(., w=30,h=12)
 network_graph$plot
dev.off()
