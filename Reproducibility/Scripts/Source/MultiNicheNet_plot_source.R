make_ligand_activity_target_plot2 = function(group_oi, order_ligands, order_targets, receiver_oi, prioritized_tbl_oi, prioritization_tables, ligand_activities_targets_DEgenes, contrast_tbl, grouping_tbl, receiver_info, ligand_target_matrix, groups_oi = NULL, plot_legend = TRUE, heights = NULL, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  if(is.null(groups_oi)){
    groups_oi = contrast_tbl %>% dplyr::pull(group) %>% unique() 
  }
  
  best_upstream_ligands = prioritized_tbl_oi$ligand %>% unique()
  
  # Ligand-Target heatmap
  active_ligand_target_links_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::select(ligand, target, ligand_target_weight, direction_regulation) %>% dplyr::rename(weight = ligand_target_weight )
  
  active_ligand_target_links_df = active_ligand_target_links_df %>% dplyr::filter(!is.na(weight))
  if(active_ligand_target_links_df$target %>% unique() %>% length() <= 2){
    cutoff = 0
  } else {
    cutoff = 0.2
  }
  
  active_ligand_target_links = nichenetr::prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = cutoff)
  
  order_ligands_ = generics::intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets_ = active_ligand_target_links_df$target %>% unique() %>% generics::intersect(rownames(active_ligand_target_links))
  
  # order_ligands = order_ligands_ %>% make.names()
  # order_targets = order_targets_ %>% make.names()
  
  rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
  
  if(!is.matrix(active_ligand_target_links[order_targets,order_ligands]) ){
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% matrix(ncol = 1)
    rownames(vis_ligand_target) = order_ligands
    colnames(vis_ligand_target) = order_targets
  } else {
    vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
  }
  
  vis_ligand_target_df = vis_ligand_target %>% data.frame() %>% tibble::rownames_to_column("ligand") %>% tidyr::gather("target","score", -ligand) %>% tibble::as_tibble() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands))  %>% dplyr::inner_join(active_ligand_target_links_df %>% distinct(target, direction_regulation)) %>% dplyr::mutate(target = factor(target, levels = order_targets))
  
  library(RColorBrewer)
  YlGn <- brewer.pal(9, 'YlGn')
  
  cutoff=0.1
  vis_ligand_target_df[vis_ligand_target_df>cutoff] = cutoff
  p_ligand_target_network = vis_ligand_target_df %>% ggplot(aes(target,ligand,fill = score)) + 
    geom_tile(color = "whitesmoke", size = 0.5) + 
    facet_grid(.~direction_regulation, scales = "free", space = "free") +
#    scale_fill_gradient2(low = "white", mid = "purple", high = "darkred", midpoint = 0.14) + 
    theme_light() +
    scale_x_discrete(position = "top") + 
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.50, "lines"),
      strip.text.x = element_text(size = 9, color = "black"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Regulatory Potential") + xlab("Predicted target genes") + ylab("Prioritized ligands")
  
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", "plum1", "orchid2","orchid4","violetred"),values = c(0, 0.05, 0.50, 0.80, 1),  limits = c(0, max(ligand_target_matrix)))
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.025, 0.075, 0.25, 0.40, 0.55, 1),  limits = c(0, max(ligand_target_matrix)))
  # custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 11, name = "PiYG") %>% .[1:5] %>% rev()),values = c(0, 0.04, 0.12, 0.30, 0.40, 0.55, 1),  limits = c(0, max(ligand_target_matrix)))
  custom_scale_fill = scale_fill_gradientn(colours = YlGn, limits = c(0, cutoff))
  
  p_ligand_target_network = p_ligand_target_network + custom_scale_fill
  
  # Ligand-Activity-Scaled -----
  ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% order_ligands & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(group %in% groups_oi) %>% dplyr::select(ligand, group, direction_regulation, activity_scaled) %>% dplyr::distinct() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands)) 
  max_activity = abs(ligand_activity_df$activity_scaled) %>% max(na.rm = TRUE)
  ligand_activity_df$activity_scaled[ligand_activity_df$activity_scaled>0.75*max_activity] <- 0.75*max_activity
  p_ligand_activity_scaled = ligand_activity_df %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(direction_regulation , ligand, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke", size = 0.5) +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(0.5, "lines"),
      strip.text.x = element_text(size = 10, color = "black"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver") + ylab("Prioritized ligands") + xlab("Scaled ligand activity")
  
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-0.75*max_activity, 0.75*max_activity))
  p_ligand_activity_scaled = p_ligand_activity_scaled + custom_scale_fill
  
  # Ligand-Activity -----
  ligand_activity_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::filter(ligand %in% order_ligands_ & receiver == receiver_oi) %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(group %in% groups_oi) %>% dplyr::select(ligand, group, direction_regulation, activity) %>% dplyr::distinct() %>% dplyr::mutate(ligand = factor(ligand, levels = order_ligands)) 
  
  p_ligand_activity = ligand_activity_df %>%
    # ggplot(aes(receiver, lr_interaction, color = activity_scaled, size = activity)) +
    # geom_point() +
    ggplot(aes(direction_regulation , ligand, fill = activity)) +
    geom_tile(color = "whitesmoke", size = 0.5) +
    facet_grid(.~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    # xlab("Ligand activities in receiver cell types\n\n") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.50, "lines"),
      panel.spacing.y = unit(0.50, "lines"),
      strip.text.x = element_text(size = 10, color = "black"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Ligand Activity\nin Receiver") + ylab("Prioritized ligands") + xlab("Ligand activity")
  custom_scale_fill = scale_fill_gradient2(low = "white", mid = "white",high = "darkorange",midpoint = 0)
  p_ligand_activity = p_ligand_activity + custom_scale_fill
  
  # Target expression
  target_regulation_df = ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::distinct(target,  direction_regulation)  %>% dplyr::rename(gene = target)
  
  # p_targets = make_DEgene_dotplot_pseudobulk_reversed(genes_oi = order_targets_, celltype_info = receiver_info, prioritization_tables = prioritization_tables, celltype_oi = receiver_oi, grouping_tbl = grouping_tbl, groups_oi = groups_oi, target_regulation_df = target_regulation_df)
  
  # -----------------------

  target_regulation_df = multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>% dplyr::ungroup() %>% dplyr::inner_join(contrast_tbl) %>% 
                       dplyr::filter(ligand %in% best_upstream_ligands & receiver == receiver_oi & group == group_oi) %>% dplyr::ungroup() %>% dplyr::distinct(target,  direction_regulation)  %>% dplyr::rename(gene = target)
  
  genes_oi = order_targets
  celltype_info = multinichenet_output$celltype_info 
  prioritization_tables = multinichenet_output$prioritization_tables 
  celltype_oi = receiver_oi 
  grouping_tbl = multinichenet_output$grouping_tbl 
  groups_oi = groups_oi 
  target_regulation_df = target_regulation_df
  
  ####  make the plot that indicates whether a celltype was sufficiently present in a sample ####
  keep_tbl = prioritization_tables$sample_prioritization_tbl %>% dplyr::distinct(sample, group, receiver, keep_receiver) %>% dplyr::rename(celltype = receiver) %>% dplyr::mutate(keep_receiver = as.logical(keep_receiver))
  
  keep_sender_receiver_values = c(1, 4)
  names(keep_sender_receiver_values) = c(FALSE, TRUE)
  
  plot_data = celltype_info$pb_df %>% dplyr::inner_join(grouping_tbl, by = c("sample")) %>% dplyr::inner_join(keep_tbl, by = c("sample","group","celltype"))
  plot_data = plot_data %>% dplyr::group_by(gene,celltype) %>% dplyr::mutate(scaled_gene_exprs = nichenetr::scaling_zscore(pb_sample)) %>% dplyr::ungroup()
  plot_data$gene = factor(plot_data$gene, levels=genes_oi)
  
  plot_data = plot_data %>% dplyr::filter(gene %in% genes_oi & celltype %in% celltype_oi) %>% dplyr::inner_join(target_regulation_df) %>% dplyr::mutate(gene = factor(gene, levels = genes_oi))
  
  if(!is.null(groups_oi)){
    plot_data = plot_data %>% dplyr::filter(group %in% groups_oi)
  }
  
  #-------------------------------------------
  RdYlBu <- rev(brewer.pal(9, 'RdYlBu'))
  plot_data$scaled_gene_exprs[-2 > plot_data$scaled_gene_exprs] = -2
  plot_data$scaled_gene_exprs[plot_data$scaled_gene_exprs > 2] = 2

  p_targets = plot_data %>%
    ggplot(aes(gene, sample, fill = scaled_gene_exprs)) + 
    geom_tile(color = "whitesmoke", size = 0.5) + 
    facet_grid(group~direction_regulation, scales = "free", space = "free") +
    # scale_fill_gradient2(low = "white", mid = "purple", high = "darkred", midpoint = 0.14) +
    scale_fill_gradientn(colors = RdYlBu,limits=c(-2,2)) +
    theme_light() +
    scale_x_discrete(position = "top") + 
    theme(
      axis.ticks = element_blank(),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0, face = "italic"),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.50, "lines"),
      strip.text.x = element_text(size = 9, color = "black"),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled pseudobulk\nexpression") + xlab("Genes") + ylab("Samples")

  # Combine the plots -----
  n_groups = ligand_activity_df$group %>% unique() %>% length()
  n_targets = ncol(vis_ligand_target)
  n_ligands = nrow(vis_ligand_target)
  n_samples = grouping_tbl %>% dplyr::filter(group %in% groups_oi) %>% dplyr::pull(sample) %>% length()
  
  legends = patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_activity_scaled)),ggpubr::as_ggplot(ggpubr::get_legend(p_ligand_target_network)), nrow = 2) %>%
    patchwork::wrap_plots(ggpubr::as_ggplot(ggpubr::get_legend(p_targets)))
  
  if(is.null(heights)){
    heights = c(n_ligands + 3, n_samples)
  }
  if(is.null(widths)){
    widths = c(n_groups*2 + 0.75, n_groups*2, n_targets)
  }
  
  if(plot_legend == FALSE){
    design <- "AaB
               ##C"
    combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_activity + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none") + xlab(""),
                                          nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
    
  } else {
    design <- "AaB
               L#C"
    
    combined_plot = patchwork::wrap_plots(A = p_ligand_activity_scaled + theme(legend.position = "none", axis.ticks = element_blank()) + theme(axis.title.x = element_text()),
                                          a = p_ligand_activity + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          B = p_ligand_target_network + theme(legend.position = "none", axis.ticks = element_blank()) + ylab(""),
                                          C = p_targets + theme(legend.position = "none") + xlab(""),
                                          L = legends, nrow = 2, design = design, widths = widths, heights = heights)
    return(list(combined_plot = combined_plot, legends = legends))
  }
}  

####################################################################################################################
####################################################################################################################

make_sample_lr_prod_activity_plots_Omnipath2 = function(prioritization_tables, prioritized_tbl_oi, widths = NULL){
  requireNamespace("dplyr")
  requireNamespace("ggplot2")
  
  sample_data = prioritization_tables$sample_prioritization_tbl %>% dplyr::filter(id %in% prioritized_tbl_oi$id) %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))   %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  sample_data = sample_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = sample_data$sender_receiver %>% unique()))
  
  group_data = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, activity, activity_scaled, direction_regulation, prioritization_score) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data_celltype_specificity = prioritization_tables$group_prioritization_tbl  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_celltype_specificity = group_data_celltype_specificity %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data_celltype_specificity$sender_receiver %>% unique()))
  
  group_data_frac_expression = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% dplyr::filter(id %in% sample_data$id) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_frac_expression = group_data_frac_expression %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  omnipath_df = prioritized_tbl_oi %>% distinct(id, curation_effort, n_references, n_resources)
  
  group_data_omnipath = prioritization_tables$group_prioritization_table_source  %>% dplyr::mutate(sender_receiver = paste(sender, receiver, sep = " --> "), lr_interaction = paste(ligand, receptor, sep = " - "))  %>% dplyr::distinct(id, sender, receiver, sender_receiver, lr_interaction) %>% dplyr::filter(id %in% sample_data$id) %>% inner_join(omnipath_df) %>%  dplyr::arrange(receiver) %>% dplyr::group_by(receiver) %>%  dplyr::arrange(sender, .by_group = TRUE)
  group_data_omnipath = group_data_omnipath %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  group_data = group_data %>% inner_join(group_data_celltype_specificity) %>% inner_join(group_data_frac_expression) 
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  rm(group_data_celltype_specificity)
  rm(group_data_frac_expression)

  keep_sender_receiver_values = c(0.25, 0.9, 1.75, 4)
  names(keep_sender_receiver_values) = levels(sample_data$keep_sender_receiver)

  p1 = sample_data %>%
    ggplot(aes(sample, lr_interaction, fill = scaled_LR_pb_prod, size = keep_sender_receiver)) +
    geom_tile(color = "whitesmoke", size = 0.5) + 
    facet_grid(sender_receiver~group, scales = "free", space = "free", switch = "y") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      # axis.title.x = element_text(face = "bold", size = 11),       axis.title.y = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.40, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x.top = element_text(size = 10, color = "black", face = "bold", angle = 0),
      strip.text.y.left = element_text(size = 9, color = "black", face = "bold", angle = 0),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled L-R\npseudobulk exprs product", size= "Sufficient presence\nof sender & receiver") + 
    scale_size_manual(values = keep_sender_receiver_values)
  max_lfc = abs(sample_data$scaled_LR_pb_prod) %>% max()
  custom_scale_fill = scale_fill_gradientn(colours = RColorBrewer::brewer.pal(n = 7, name = "RdBu") %>% rev(),values = c(0, 0.350, 0.4850, 0.5, 0.5150, 0.65, 1),  limits = c(-1*max_lfc, max_lfc))
  
  p1 = p1 + custom_scale_fill
  

  p2 = group_data %>% dplyr::filter(., activity_scaled > 0.5) %>%
    ggplot(aes(direction_regulation , lr_interaction, fill = activity_scaled)) +
    geom_tile(color = "whitesmoke") +
    facet_grid(sender_receiver~group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_text(face = "bold.italic", size = 9),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(fill = "Scaled Ligand\nActivity in Receiver")
  max_activity = abs(group_data$activity_scaled) %>% max(na.rm = TRUE)
  custom_scale_fill = scale_fill_gradientn(colours = c("white", RColorBrewer::brewer.pal(n = 7, name = "PuRd") %>% .[-7]),values = c(0, 0.51, 0.575, 0.625, 0.675, 0.725, 1),  limits = c(-1*max_activity, max_activity))
  
  p2 = p2 + custom_scale_fill
  
  # add the plot visualizing cell-type specificity
  # cs_data = group_data %>% filter(group %in% prioritized_tbl_oi$group) %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data = group_data %>% distinct(sender_receiver, lr_interaction, group, scaled_pb_ligand, scaled_pb_receptor) %>% tidyr::gather(LR, celltype_specificity, scaled_pb_ligand:scaled_pb_receptor)
  cs_data$LR[cs_data$LR == "scaled_pb_ligand"] = "ligand"
  cs_data$LR[cs_data$LR == "scaled_pb_receptor"] = "receptor"
  frac_data = group_data %>% distinct(sender_receiver, lr_interaction, group, fraction_ligand_group, fraction_receptor_group) %>% tidyr::gather(LR, fraction_expression, fraction_ligand_group:fraction_receptor_group)
  frac_data$LR[frac_data$LR == "fraction_ligand_group"] = "ligand"
  frac_data$LR[frac_data$LR == "fraction_receptor_group"] = "receptor"
  
  cs_data = cs_data %>% inner_join(frac_data)
  
  p_cs = cs_data %>% 
    ggplot(aes(LR , lr_interaction, color = celltype_specificity, size = fraction_expression)) +
    geom_point() +
    facet_grid(sender_receiver ~ group, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    viridis::scale_color_viridis() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    ) + labs(color = "Scaled celltype specificity") + labs(size = "Fraction of expression")
  
  # add the plot visualizing Omnipath DB LR scores
  group_data = group_data %>% inner_join(group_data_omnipath) 
  group_data = group_data %>% dplyr::mutate(sender_receiver = factor(sender_receiver, levels = group_data$sender_receiver %>% unique()))
  
  omnipath_data = group_data %>% distinct(sender_receiver, lr_interaction, curation_effort, n_references, n_resources) %>% tidyr::gather(omnipath_score_type, omnipath_score, curation_effort:n_resources)
  omnipath_data$omnipath_score_type[omnipath_data$omnipath_score_type == "curation_effort"] = "curation effort"
  omnipath_data$omnipath_score_type[omnipath_data$omnipath_score_type == "n_references"] = "nr. references"
  omnipath_data$omnipath_score_type[omnipath_data$omnipath_score_type == "n_resources"] = "nr. resources"
  omnipath_data = omnipath_data %>% mutate(omnipath = "Omnipath")
  p_omnipath = omnipath_data %>% 
    ggplot(aes(omnipath_score_type , lr_interaction, size = omnipath_score)) +
    geom_point(color = "grey40") +
    facet_grid(sender_receiver ~ omnipath, scales = "free", space = "free") +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 9,  angle = 90,hjust = 0),
      strip.text.x.top = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(0.20, "lines"),
      panel.spacing.y = unit(0.25, "lines"),
      strip.text.x = element_text(size = 10, color = "black", face = "bold"),
      strip.text.y = element_blank(),
      strip.background = element_rect(color="darkgrey", fill="whitesmoke", size=1.5, linetype="solid")
    )  + labs(size = "Omnipath DB score")
  
  
  if(!is.null(widths)){
    p = patchwork::wrap_plots(
      p1,p2,p_cs, p_omnipath,
      nrow = 1,guides = "collect",
      widths = widths
    )
  } else {
    p = patchwork::wrap_plots(
      p1,p2,p_cs, p_omnipath,
      nrow = 1,guides = "collect",
      widths = c(sample_data$sample %>% unique() %>% length(), 2*(sample_data$group %>% unique() %>% length()),2*(sample_data$group %>% unique() %>% length()),3)
    )
  }
  
  return(p)
  
}