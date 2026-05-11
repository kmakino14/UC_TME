#####################################################
# TREKKER ArchR NFR
#####################################################

setwd(path_to_wd)
source("Reproducibility/Scripts/Source/my.source.R")

suppressMessages(library(ArchR))
suppressMessages(library(parallel))
suppressMessages(library(dplyr))
set.seed(1234)

SAMPLE_TMP = 'UC_TREKKER'

######################################################
# Create Arrow files
######################################################

fragment_path  = paste0("Reproducibility/Data/ATAC_fragments/atac_fragments_TREKKER.tsv.gz")   #### cellranger-atac aggr output after subsetting from ARC

dir_list = c(paste0("Reproducibility/Results/ArchR/TREKKER/output/log"), 
             paste0("Reproducibility/Results/ArchR/TREKKER/output/QC")
             )
for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}

addArchRThreads(threads = 1) 
addArchRGenome("hg38")

ArrowFiles = createArrowFiles(
  inputFiles      = fragment_path,
  sampleNames     = SAMPLE_TMP,
  minTSS       = 1, #Dont set this too high because you can always increase later
  minFrags     = 500, 
  addTileMat      = TRUE,
  addGeneScoreMat = FALSE,
  force = TRUE,
  QCDir   = "Reproducibility/Results/ArchR/TREKKER/output/QC", 
  logFile = "Reproducibility/Results/ArchR/TREKKER/output/log/UC_DOGMA_createArrowFiles.log"
)

proj = ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = paste0("Reproducibility/Results/ArchR/TREKKER/output"),
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)

saveArchRProject(ArchRProj = proj, 
                 outputDirectory = paste0("Reproducibility/Results/ArchR/TREKKER/output"),
                 overwrite = FALSE, load = TRUE, dropCells = FALSE,
                 logFile = createLogFile("saveArchRProject"),
                 threads = 1
                 )

######################################################
# ArchR project construction by lineage
######################################################

ArchR_preprocess = function(tmp_lineage){
    proj = loadArchRProject(path = paste0("Reproducibility/Results/ArchR/TREKKER/output"),
                            force = FALSE, showLogo = FALSE)

    metadata = fread_n("Reproducibility/Data/TREKKER/UC_TREKKER_metadata.txt")
    df = dplyr::filter(metadata, celltype %in% 'Epithelial')
    
    ######################################################
    
    dir_list = c(paste0("Reproducibility/Results/ArchR/",tmp_lineage,"/output/log"),
                 paste0("Reproducibility/Results/ArchR/",tmp_lineage,"/output/QC"),
                 paste0("Reproducibility/Results/ArchR/",tmp_lineage,"/export")
                 )
    for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}
        
    ## cell barcode after QC filter by Seurat
    target_cell = intersect(paste0(SAMPLE_TMP,"#",rownames(df)),rownames(getCellColData(proj)))
    proj_subset = subsetArchRProject(
      ArchRProj = proj,
      cells = target_cell,
      outputDirectory = paste0("Reproducibility/Results/ArchR/",tmp_lineage,"/output"),
      dropCells = TRUE,
      logFile = NULL,
      threads = getArchRThreads(),
      force = TRUE
    )

    ######################################################
    # Preprocesing
    # SVD, Clustering, UMAP
    proj_subset <- addIterativeLSI(ArchRProj = proj_subset, useMatrix = "TileMatrix", name = "IterativeLSI", force=TRUE)
    
    # Gene scores with selected features
    # Artificial black list to exclude all non variable features
    chrs <- getChromSizes(proj_subset)
    var_features <- proj_subset@reducedDims[["IterativeLSI"]]$LSIFeatures
    var_features_gr <- GRanges(var_features$seqnames, IRanges(var_features$start, var_features$start + 500))
    blacklist <- setdiff(chrs, var_features_gr)
    proj_subset <- addGeneScoreMatrix(proj_subset, matrixName='GeneScoreMatrix', force=TRUE, blacklist=blacklist)

    # Peaks using NFR fragments
    # proj_subset <- addClusters(input = proj_subset, reducedDims = "IterativeLSI")
    target_cell = intersect(paste0(SAMPLE_TMP,"#",rownames(df)),rownames(getCellColData(proj_subset)))
    df$CB = paste0(SAMPLE_TMP,"#",rownames(df))
    Single_cell_list = left_join(data.frame(CB=target_cell), df ,by="CB")
    proj_subset$celltype = Single_cell_list$clone
    
    proj_subset <- addGroupCoverages(proj_subset, maxFragmentLength=147, groupBy = "celltype")
    proj_subset <- addReproduciblePeakSet(proj_subset, groupBy = "celltype")
    
    # Counts
    proj_subset <- addPeakMatrix(proj_subset, maxFragmentLength=147, ceiling=10^9)
    
    # Save 
    proj_subset <- saveArchRProject(ArchRProj = proj_subset)
    
    ######################################################
    # Export
    write.csv(getReducedDims(proj_subset), paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/svd.csv'), quote=FALSE)
    write.csv(getCellColData(proj_subset), paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/cell_metadata.csv'), quote=FALSE)
    
    # Gene scores
    gene.scores <- getMatrixFromProject(proj_subset)
    scores <- assays(gene.scores)[['GeneScoreMatrix']]
    scores <- as.matrix(scores)
    rownames(scores) <- rowData(gene.scores)$name
    write.csv(scores, paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/gene_scores.csv'), quote=FALSE)
    
    # Peak counts
    peaks <- getPeakSet(proj_subset)
    peak.counts <- getMatrixFromProject(proj_subset, 'PeakMatrix')
    
    # Reorder peaks 
    # Chromosome order
    chr_order <- sort(seqlevels(peaks))
    reordered_features <- list()
    for(chr in chr_order)
        reordered_features[[chr]] = peaks[seqnames(peaks) == chr]
    reordered_features <- Reduce("c", reordered_features)    
    
    # Export counts
    dir_list = c(paste0("Reproducibility/Results/ArchR/",tmp_lineage,"/export/peak_counts"))
    for(dir_name_tmp in dir_list ){dir.create(dir_name_tmp,   showWarnings = FALSE, recursive = TRUE)}
    
    counts <- assays(peak.counts)[['PeakMatrix']]
    writeMM(counts, paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/peak_counts/counts.mtx'))
    write.csv(colnames(peak.counts), paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/peak_counts/cells.csv'), quote=FALSE)
    names(reordered_features) <- sprintf("Peak%d", 1:length(reordered_features))
    write.csv(as.data.frame(reordered_features), paste0('Reproducibility/Results/ArchR/',tmp_lineage,'/export/peak_counts/peaks.csv'), quote=FALSE)
}

ArchR_preprocess(tmp_lineage = 'TREKKER_Epithelial')