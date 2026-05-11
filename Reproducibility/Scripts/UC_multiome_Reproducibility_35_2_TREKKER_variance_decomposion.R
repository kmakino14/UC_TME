#####################################################
# Variance decomposition
#####################################################

setwd(path_to_wd)
source("~/my.source.R")
options(stringsAsFactors=F)
source("Reproducibility/Scripts/Source/my.source.R")
source("Reproducibility/Scripts/Source/Seurat_source.R")
options(stringsAsFactors=F)

suppressMessages(library(Signac))
suppressMessages(library(Seurat))
suppressMessages(library(SingleCellExperiment))
suppressMessages(library(cicero))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(VISION))

suppressMessages(library(GenomeInfoDb))
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg38))

suppressMessages(library(ggplot2))
suppressMessages(library(ggthemes))
suppressMessages(library(ggrepel))
suppressMessages(library(ggbeeswarm))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(circlize))

suppressMessages(library(viridis))
suppressMessages(library(RColorBrewer))
suppressMessages(library(ggsci))
suppressMessages(library(colorspace))
suppressMessages(library(BuenColors))

suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(escape))
suppressMessages(library(forestmodel))

suppressMessages(library(dplyr))
suppressMessages(library(jsonlite))

set.seed(1234)

####################################################################################################
## P02
####################################################################################################
#==========
# Genes
#==========

df_anova = fread_n("Reproducibility/Results/TREKKER/scVIVA/P02/model/hypoxia_niche_anova_df_15.txt")
df_anova <- df_anova[complete.cases(df_anova), ]

xvals <- df_anova$Signed_Clone_Effect
yvals <- df_anova$Signed_Env_Effect
lims <- range(c(xvals, yvals), na.rm = TRUE)

hypoxia_genes <- c(
  "ADM", "AHNAK", "ALDOA", "ANGPTL4", "BNIP3", "BNIP3L", "C4orf3", "BLNK", 
  "CA9", "CASP14", "DDIT4", "DHRS3", "ENO2",  "ERO1L", "EPAS1","EGLN3", 
  "FAM162A", "FUT11", "HILPDA", "HK2", "IGFBP3", "LDHA", "LINC01133", "MGP",
  "NDRG1", "NDUFA4L2", "NRN1", "P4HA1", "PDK1", "PFKFB3", "PFKP", "PGF",
  "PGK1", "PLIN2", "PLOD2", "PTHLH", "RGCC", "SLC16A1", "SLC16A3", "SLC2A1",
  "SLC2A3", "SLC6A8"
)

pEMT_genes <- c(
"COL17A1"   , "NDRG1"   ,   "SFN"       , "NRG1"      , "RAF1"     , 
"KRT14"     , "KRT5"    ,   "LAMC2"     , "VGLL4"     , "IQSEC1"   , 
"COL5A1"    , "CAV1"    ,   "TGFBI"     , "F3"        , "BNC1"     , 
"MT2A"      , "PTHLH"   ,   "VSTM2L"    , "HEPHL1"    ,             
"SLC1A6"    , "CALML5"  ,   "HLA-G"     , "CDKN2B"    , "MSLN"     , 
"FST"       , "CAND2"   ,   "FKBP5"     , "FLRT2"     , "STXBP6"   , 
"PPP2R2C"   , "MUC16"   ,   "TUBB6"     , "ITGA5"     , "SERPINB2" , 
"COLEC12"   , "TMCC3"   ,   "CD109"     , "GALNT18"   , "CALB2"    , 
"CASP14"    , "CITED4"  ,   "VEGFC"     , "CXCL14"    , "CSF3R"    , 
"FGFBP1"    , "SLAMF7"  ,   "MYO16-AS1" , "SLC47A2"   ,  "HAS2"    , 
"GJA5"      , "ALPP"    ,   "SERPINB7"  , "MEI1"      ,  "KL"      ,         
"CSF2"      , "RXFP1"   ,   "KLK5"      , "MS4A15"    , "TM4SF19"  , 
"IGFL2-AS1" , "KLK7"    ,   "SERPINB9"  , "ZNF114"    , "FAP"      , 
"TRIML2"    , "SLCO2B1" ,   "RGS20"     , "KLK6"      , "HAS2-AS1" , 
"TLL1"      , "KRT81"   ,   "APOB"      , "KRT6C"     , "HKDC1"    , 
"MFAP5"     , "JCAD"    ,   "ECM1"      , "CST6"      , "CD274"    , 
"INSL4"     , "IGFL1"   ,   "SLC34A2"   , "ANPEP"     , "PSG11"    , 
"ALPG"      , "LRP2"    ,   "TRIB3" 
)

urothelium <- c("UPK1A", "UPK3A", "KRT20", "PSCA", "SNX31", "PPARG", 
                "KRT5", "KRT14", 'KRT17', "TP63", "CD44", "EGFR",
                "GATA3", "KRT13", "EPCAM")

hypoxia_genes <- intersect(hypoxia_genes, df_anova$Gene)
pEMT_genes <- intersect(pEMT_genes, df_anova$Gene)
urothelium <- intersect(urothelium, df_anova$Gene)

Create_scatter_plot <- function(df, highlight_genes, limits, target_genes, color){
  df_highlight <- df %>% dplyr::filter(Gene %in% highlight_genes)
  df_others    <- df %>% dplyr::filter(!Gene %in% highlight_genes)

  df_highlight <- df_highlight %>%
  mutate(label = ifelse(Gene %in% target_genes, as.character(Gene), NA_character_))

  ggplot() +
   geom_point(
     data = df_others,
     aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
     color = "grey75", 
     alpha = 0.3, 
     size = 0.5, 
     shape = 16
   ) +
   # y = 0
   geom_hline(
     yintercept = 0,
     linetype = "dotted", color = "black", linewidth = 0.5
   ) +
   geom_vline(
     xintercept = 0,
     linetype = "dotted", color = "black", linewidth = 0.5
   ) +
   geom_point(
     data = df_highlight,
     aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
     color = color, 
     alpha = 0.9, 
     size = 2.0, 
     shape = 16
   ) +
   geom_text_repel(
     data = df_highlight,
     aes(x = Signed_Clone_Effect, y = Signed_Env_Effect, label = label),
     size = 3,              
     box.padding = 0.5,     
     point.padding = 0.3,  
     max.overlaps = Inf,   
     color = "black",      
     fontface = "bold",    
     min.segment.length = 0 
   ) +
   coord_fixed(xlim = limits, ylim = limits) +
   theme_classic(base_size = 14) +
   labs(
     x = "Signed Clone Effect",
     y = "Signed Environmental Effect"
   )
}

target_genes = c('CA9','SLC2A1','LDHA', 'PGK1','PDK1')
p1 = Create_scatter_plot(df_anova, hypoxia_genes, lims, target_genes, "#de2d26")

target_genes = c('NDRG1','VEGFC','CDKN2B','COL17A1')
p2 = Create_scatter_plot(df_anova, pEMT_genes, lims, target_genes, "#6A3D9A")

target_genes = c('GATA3','PPARG','UPK1A','KRT20', 'KRT17')
p3 = Create_scatter_plot(df_anova, urothelium, lims, target_genes,  "#008B8B")

paste0("Reproducibility/Results/Plots/Slide-tags/Figure3L_genes.pdf") %>% pdf(w=15, h=4)
 plot(p1|p2|p3)
dev.off()

#==========
# TFs
#==========

df_anova = fread_n("Reproducibility/Results/TREKKER/scVIVA/P02/model/hypoxia_niche_anova_df_15_TF.txt")
df_anova <- df_anova[complete.cases(df_anova), ]

xvals <- df_anova$Signed_Clone_Effect
yvals <- df_anova$Signed_Env_Effect
lims <- range(c(xvals, yvals), na.rm = TRUE)

Create_scatter_plot <- function(df, highlight_genes, limits, target_genes, color1, color2){
  df_highlight1  <- df %>% dplyr::filter(Gene %in% highlight_genes)
  df_highlight2 <- df %>% dplyr::filter(Gene %in% target_genes)
  df_others     <- df %>% dplyr::filter(!Gene %in% c(highlight_genes, target_genes))

  df_highlight <- rbind(df_highlight1, df_highlight2) %>%
    dplyr::mutate(label = as.character(Gene))

  ggplot() +
    geom_point(
      data = df_others,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
      color = "grey75",
      alpha = 0.3,
      size = 1.6,
      shape = 16
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dotted", color = "black", linewidth = 0.5
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dotted", color = "black", linewidth = 0.5
    ) +
    geom_point(
      data = df_highlight1,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
      color = color1,
      alpha = 0.9,
      size = 2.0,
      shape = 16
    ) +
    geom_point(
      data = df_highlight2,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
      color = color2,
      alpha = 0.95,
      size = 2.2,
      shape = 16
    ) +
    geom_text_repel(
      data = df_highlight,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect, label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = Inf,
      color = "black",
      fontface = "bold",
      min.segment.length = 0
    ) +
    coord_fixed(xlim = limits, ylim = limits) +
    theme_classic(base_size = 14) +
    labs(
      x = "Signed Clone Effect",
      y = "Signed Environmental Effect"
    )
}

lum_TF = c('GATA3','FOXA1')
env_TF = c('SMAD3','HIF1A', 'MYCN')
p1 = Create_scatter_plot(df_anova, env_TF, lims, lum_TF, "#de2d26", "#008B8B")  #  

paste0("Reproducibility/Results/Plots/Slide-tags/Figure3L_TF.pdf") %>% pdf(w=4, h=4)
 plot(p1)
dev.off()

####################################################################################################
## P06
####################################################################################################
#==========
# Genes
#==========

df_anova = fread_n("Reproducibility/Results/TREKKER/scVIVA/P06/model/stress_niche_anova_df_19.txt")
df_anova <- df_anova[complete.cases(df_anova), ]

stress_genes <- c(
  "ABHD3", "AC016629.8", "ADAMTS1", "ANKRD28", "ANKRD37", "ARC", "AREG",
  "ARID5B", "ARL5B", "ATF3", "B3GNT5", "BAG3", "BAMBI", "BBC3", "BHLHE40",
  "BIRC3", "BRD2", "BTG1", "BTG2", "C10orf10", "C11orf96", "C16orf98",
  "C1orf63", "C5orf45", "C8orf4", "CA12", "CCL20", "CCNL1", "CCP110",
  "CD55", "CDKN1A", "CDKN2AIP", "CEBPB", "CEBPD", "CITED2", "CLDN3",
  "CLDN4", "CLK1", "CSRNP1", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL8",
  "CYP1B1", "CYR61", "DAPP1", "DEPP1", "DNAJA1", "DNAJA4", "DNAJB1",
  "DNAJB4", "DSC3", "DUSP1", "DUSP10", "DUSP14", "DUSP2", "DUSP5", "EDN1",
  "EFNA1", "EGR1", "EGR2", "EGR3", "EHF", "EIF4A2", "ELF3", "EMP1",
  "EPHA2", "EREG", "ERRFI1", "ETS2", "FADS3", "FAM133B", "FAM46A",
  "FAM53C", "FGFR3", "FHL2", "FKBP4", "FOS", "FOSB", "FOSL1", "FOSL2",
  "FRMD4B", "FUS", "GADD45A", "GADD45B", "GADD45G", "GDF15", "GEM",
  "GOLGB1", "HBEGF", "HES1", "HEXIM1", "HNRNPU-AS1", "HSP90AA1", "HSPA1A",
  "HSPA1B", "HSPA5", "HSPA6", "HSPA8", "HSPB1", "HSPH1", "ID1", "ID2",
  "IER2", "IER3", "IER5", "IFRD1", "IL8", "IRX2", "ITPKC", "JAG1", "JUN",
  "JUNB", "JUND", "KCNQ1OT1", "KDM6B", "KIAA1683", "KLF10", "KLF2", "KLF4",
  "KLF5", "KLF6", "LDLR", "LIF", "LIMA1", "LINC00910", "LINC00936", "LMNA",
  "MAFB", "MAFF", "MALAT1", "MAP3K8", "MCL1", "METTL15", "MIDN", "MIR222HG",
  "MIR24-2", "MIR3064", "MKNK2", "MXD1", "MYADM", "MYC", "MYLIP", "NAMPT",
  "NCOA7", "NEAT1", "NEDD9", "NFKBIA", "NFKBID", "NFKBIZ", "NOP58", "NPAS4",
  "NR4A1", "NR4A2", "NRARP", "NSUN6", "OVOL1", "PDK4", "PDLIM3", "PGM2L1",
  "PHLDA1", "PHLDA2", "PIM1", "PIM3", "PLAUR", "PLK2", "PLK3", "PMAIP1",
  "PNRC1", "POLG2", "PPP1R10", "PPP1R15A", "PRKCA", "PRPF40A", "PTGS2",
  "PTRHD1", "RAB11FIP1", "RAB3IP", "RALGDS", "RASD1", "RDH10", "REL", "RGS2",
  "RHOB", "RND1", "RND3", "RP1-313I6.12", "RP11-182L21.6", "RP5-821D11.7",
  "RP6-99M1.2", "RRAD", "SEMA4A", "SERTAD1", "SFN", "SIK1", "SLC20A1",
  "SLC38A2", "SNHG12", "SOCS1", "SOCS3", "SOX4", "SOX9", "SQSTM1", "TACSTD2",
  "TAF1D", "TCF7", "TCIM", "TFRC", "TGIF1", "TIPARP", "TLE4", "TM4SF1",
  "TNFAIP2", "TNFAIP3", "TRIB1", "TSC22D1", "TSC22D3", "UBC", "UGDH",
  "VPS37B", "WEE1", "YBX3", "YME1L1", "YOD1", "ZBTB43", "ZC3H12A", "ZFAND2A",
  "ZFP36", "ZFP36L1", "ZFP36L2", "ZNF165", "ZNF430"
)

stress_genes <- intersect(stress_genes, df_anova$Gene)

df_bg   <- subset(df_anova, !(Gene %in% c(stress_genes)))
df_str  <- subset(df_anova, Gene %in% stress_genes)

xvals <- df_anova$Signed_Clone_Effect
yvals <- df_anova$Signed_Env_Effect
lims <- range(c(xvals, yvals), na.rm = TRUE)

create_base_plot_with_density <- function(bg_data, limits, label_genes = NULL, max_density = 20) {
  dens <- MASS::kde2d(bg_data$Signed_Clone_Effect, bg_data$Signed_Env_Effect, n = 100)
  ix <- findInterval(bg_data$Signed_Clone_Effect, dens$x)
  iy <- findInterval(bg_data$Signed_Env_Effect, dens$y)
  ii <- cbind(ix, iy)
  bg_data$density_val <- dens$z[ii]
  
  p <- ggplot(bg_data, aes(x = Signed_Clone_Effect, y = Signed_Env_Effect)) +
    geom_point_rast(
      aes(color = density_val),
      size = 0.5,
      alpha = 0.6,
      shape = 16,
      raster.dpi = 300
    ) +
    scale_color_gradientn(
      colors = jdb_palette("horizon"), 
      limits = c(0, max_density),
      oob = scales::squish,
      name = "Density"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +
    coord_fixed(xlim = limits, ylim = limits) +
    theme_classic(base_size = 14) +
    labs(title = "Scatter plot all genes",
         x = "Signed Clone Effect", y = "Signed Environmental Effect")

  if (!is.null(label_genes) && length(label_genes) > 0) {
    label_data <- bg_data[bg_data$Gene %in% label_genes, ]
    
    p <- p +
      geom_point(
        data = label_data,
        aes(color = density_val), 
        size = 1.5,
        shape = 16
      ) +
      geom_text_repel(
        data = label_data,
        aes(label = Gene),
        size = 4,
        color = "black",
        fontface = "italic",
        box.padding = 0.5,
        point.padding = 0.3,
        segment.color = "grey50",
        max.overlaps = Inf
      )
  }
  
  return(p)
}


create_highlight_plot <- function(bg_data, highlight_genes, limits, target_genes, max_density = 20) {
  dens <- MASS::kde2d(bg_data$Signed_Clone_Effect, bg_data$Signed_Env_Effect, n = 100)
  ix <- findInterval(bg_data$Signed_Clone_Effect, dens$x)
  iy <- findInterval(bg_data$Signed_Env_Effect, dens$y)
  ii <- cbind(ix, iy)
  bg_data$density_val <- dens$z[ii]

  df_highlight <- bg_data %>% dplyr::filter(Gene %in% highlight_genes)
  df_others    <- bg_data %>% dplyr::filter(!Gene %in% highlight_genes)
  df_highlight <- df_highlight %>%
  mutate(label = ifelse(Gene %in% target_genes, as.character(Gene), NA_character_))

  ggplot() +
    geom_point_rast(
      data = df_others,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect),
      color = "grey90", 
      size = 0.5,
      alpha = 0.3,
      shape = 16,
      raster.dpi = 300
    ) +  
    geom_point(
      data = df_highlight,
      aes(x = Signed_Clone_Effect, y = Signed_Env_Effect, color = density_val),
      size = 1.5,  
      alpha = 0.8,
      shape = 16
    ) +
    geom_text_repel(
    data = df_highlight,
    aes(x = Signed_Clone_Effect, y = Signed_Env_Effect, label = label),
    size = 3,             
    box.padding = 0.5,    
    point.padding = 0.3,  
    max.overlaps = Inf,    
    color = "black",      
    fontface = "bold",     
    min.segment.length = 0 
    ) +
    scale_color_gradientn(
      colors = jdb_palette("horizon"), 
      limits = c(0, max_density),
      oob = scales::squish,
      name = "Density"
    ) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey40", linewidth = 0.4) +
    coord_fixed(xlim = limits, ylim = limits) +
    theme_classic(base_size = 14) +
    labs(
      title = paste("Scatterplot stress_genes (n =", nrow(df_highlight), ")"),
      x = "Signed Clone Effect", 
      y = "Signed Environmental Effect"
    )
}

p <- create_base_plot_with_density(
  bg_data = df_anova, 
  limits = lims, 
  max_density = 5
)

target_genes <- c("ATF3", "HSP90AA1", "FOS", "DUSP1", "KLF4")
p_stress_only <- create_highlight_plot(df_anova, stress_genes , lims, target_genes, max_density = 5)

paste0("Reproducibility/Results/Plots/Slide-tags/FigureS5J_density.pdf") %>% pdf(w=10, h=4)
 plot(p|p_stress_only)
dev.off()

#==========
# TFs
#==========

df_anova = fread_n("Reproducibility/Results/TREKKER/scVIVA/P06/model/stress_niche_anova_df_19_TF.txt")
df_anova <- df_anova[complete.cases(df_anova), ]

# 軸範囲を共通化
xvals <- df_anova$Signed_Clone_Effect
yvals <- df_anova$Signed_Env_Effect
lims <- range(c(xvals, yvals), na.rm = TRUE)

lum_TF = c('GATA3','FOXA1')
env_TF = c('KLF4','FOS','JUN','NFE2L2', 'ATF4')
p1 = Create_scatter_plot(df_anova, env_TF, lims, lum_TF , "#de2d26", "#008B8B")  #  

paste0("Reproducibility/Results/Plots/Slide-tags/FigureS5J_TF.pdf") %>% pdf(w=4, h=4)
 plot(p1)
dev.off()