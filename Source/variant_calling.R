library(SummarizedExperiment)
library(Matrix)
library(dplyr)
library(data.table)
"%ni%" <- Negate("%in%")

call_mutations_mgatk <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10){
  
  # Determinie key coverage statistics every which way
  cov <- assays(SE)[["coverage"]]
  ref_allele <- toupper(as.character(rowRanges(SE)$refAllele))
  
  # Process mutation for one alternate letter
  process_letter <- function(letter){
    print(letter)
    boo <- ref_allele != letter & ref_allele != "N"
    pos <- start(rowRanges(SE))
    variant_name <- paste0(as.character(pos), ref_allele, ">", letter)[boo]
    nucleotide <- paste0(ref_allele, ">", letter)[boo]
    position_filt <- pos[boo]
    
    # Single cell functions
    getMutMatrix <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_fw  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_fw")]]) / cov_fw)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    getMutMatrix_rev  <- function(letter){
      mat <- ((assays(SE)[[paste0(letter, "_counts_rev")]]) / cov_rev)[boo,]
      rownames(mat) <- variant_name
      mat <- as(mat, "dgCMatrix")
      return(mat)
    }
    
    # Bulk functions
    getBulk <- function(letter){
      vec <- (Matrix::rowSums(assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / Matrix::rowSums(cov))[boo]
      return(vec)
    }
    rowVars <- function(x, ...) {
      Matrix::rowSums((x - Matrix::rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
    }
    
    update_missing_w_zero <- function(vec){
      ifelse(is.na(vec)  | is.nan(vec), 0, vec)
    }
    # Set up correlation per non-zero mutation based on the strands
    dt <- merge(data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_fw")]][boo,])), 
                data.table(Matrix::summary(assays(SE)[[paste0(letter, "_counts_rev")]][boo,])), 
                by.x = c("i", "j"), by.y = c("i", "j"), 
                all = TRUE)[x.x >0 | x.y >0]
    dt$x.x <- update_missing_w_zero(dt$x.x)
    dt$x.y <- update_missing_w_zero(dt$x.y)
    
    dt2 <- data.table(variant = variant_name[dt[[1]]],
                      cell_idx = dt[[2]], 
                      forward = dt[[3]],
                      reverse = dt[[4]])
    rm(dt)
    cor_dt <- dt2[, .(cor = cor(c(forward), c(reverse), method = "pearson", use = "pairwise.complete")), by = list(variant)]
    
    # Put in vector for convenience
    cor_vec_val <- cor_dt$cor
    names(cor_vec_val) <- as.character(cor_dt$variant )
    
    # Compute the single-cell data
    mat <- getMutMatrix(letter)
    mmat <- sparseMatrix(
      i = c(summary(mat)$i,dim(mat)[1]),
      j = c(summary(mat)$j,dim(mat)[2]),
      x = c(update_missing_w_zero(summary(mat)$x), 0)
    )
    
    # Compute bulk statistics
    mean = update_missing_w_zero(getBulk(letter))
    
    # Stablize variances by replacing low coverage cells with mean
    if(stabilize_variance){
      
      # Get indices of cell/variants where the coverage is low and pull the mean for that variant
      idx_mat <- which(data.matrix(cov[boo,] < low_coverage_threshold), arr.ind = TRUE)
      idx_mat_mean <- mean[idx_mat[,1]]
      
      # Now, make sparse matrices for quick conversion
      ones <- 1 - sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = 1
      )
      
      means_mat <- sparseMatrix(
        i = c(idx_mat[,1], dim(mmat)[1]),
        j = c(idx_mat[,2], dim(mmat)[2]),
        x = c(idx_mat_mean, 0)
      )
      
      mmat2 <- mmat * ones + means_mat
      variance = rowVars(mmat2)
      rm(mmat2); rm(ones); rm(means_mat); rm(idx_mat); rm(idx_mat_mean)
      
    } else {
      variance = rowVars(mmat)
    }
    
    detected <- (assays(SE)[[paste0(letter, "_counts_fw")]][boo,] >= 2) + (assays(SE)[[paste0(letter, "_counts_rev")]][boo,] >=2 )
    
    # Compute per-mutation summary statistics
    var_summary_df <- data.frame(
      position = position_filt,
      nucleotide = nucleotide, 
      variant = variant_name,
      vmr = variance/(mean + 0.00000000001),
      mean = round(mean,7),
      variance = round(variance,7),
      n_cells_conf_detected = Matrix::rowSums(detected == 2),
      n_cells_over_5 = Matrix::rowSums(mmat >= 0.05), 
      n_cells_over_10 = Matrix::rowSums(mmat >= 0.10),
      n_cells_over_20 = Matrix::rowSums(mmat >= 0.20),
      strand_correlation = cor_vec_val[variant_name],
      mean_coverage = Matrix::rowMeans(cov)[boo], 
      stringsAsFactors = FALSE, row.names = variant_name
    )
    se_new <- SummarizedExperiment(
      rowData = var_summary_df, 
      colData = colData(SE), 
      assays = list(allele_frequency = mmat, coverage = cov[boo,])
    )
    return(se_new)
  }
  
  return(SummarizedExperiment::rbind(process_letter("A"), process_letter("C"), process_letter("G"), process_letter("T")))
  
}

# Peter van Galen, 210215
# General functions for analyses in the MAESTER project

# General
message("cutf()")
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Function that computes all heteroplasmic variants from MAEGATK output (from Caleb Lareau). 
# Rows represents a position along the mitochondrial genome and the three possible disagreements with the reference
# (except 3107 has four possible disagreements because the reference is N)
#  message("computeAFMutMatrix()")
# computeAFMutMatrix <- function(SE){
#     cov <- assays(SE)[["coverage"]]+ 0.000001
#     ref_allele <- as.character(rowRanges(SE)$refAllele)
#     
#     getMutMatrix <- function(letter){
#       mat <- (assays(SE)[[paste0(letter, "_counts_fw")]] + assays(SE)[[paste0(letter, "_counts_rev")]]) / cov
#       rownames(mat) <- paste0(as.character(1:dim(mat)[1]), "_", toupper(ref_allele), ">", letter)
#       return(mat[toupper(ref_allele) != letter,])
#     }
#     
#     rbind(getMutMatrix("A"), getMutMatrix("C"), getMutMatrix("G"), getMutMatrix("T"))
#   }

computeAFMutMatrix <- function(SE, chromosome_prefix = "chrM"){
  cov <- assays(SE)[["coverage"]] + 0.000001
  ref_allele <- as.character(rowRanges(SE)$refAllele)

  getMutMatrix <- function(letter){
    names_rows <- paste0(chromosome_prefix, "_", 1:nrow(cov), "_", toupper(ref_allele), "_", letter)
    names_rows <- names_rows[toupper(ref_allele) != letter]
    mat_fow <- assays(SE)[[paste0(letter, "_counts_fw")]]
    mat_rev <- assays(SE)[[paste0(letter, "_counts_rev")]]
    mat <- mat_fow + mat_rev
    mat <- mat[toupper(ref_allele) != letter,]
    cov_use <- cov[toupper(ref_allele) != letter,]
    mat <- mat / cov_use
    gc()
    mat[is.na(mat)] <- 0
    # We can get AF values greater than 1, which is due to uninformative reads.
    # See: https://gatk.broadinstitute.org/hc/en-us/articles/360035532252-Allele-Depth-AD-is-lower-than-expected
    # and https://github.com/caleblareau/mgatk/issues/1
    # We simply set these values to 1, since that is the actual information we have in this case.
    # This issue can be solved on the MAEGATK/GATK side.
    mat[mat > 1] <- 1
    rownames(mat) <- names_rows
    #mat <- as(mat, "dgCMatrix")
    mat <- as(mat, "CsparseMatrix")
    return(mat)
  }
   A_matrix <- getMutMatrix("A")
   #A_matrix <- as.matrix(A_matrix)
   gc()
   C_matrix <- getMutMatrix("C")
   #C_matrix <- as.matrix(C_matrix)
   gc()
   G_matrix <- getMutMatrix("G")
   #G_matrix <- as.matrix(G_matrix)
   gc()
   T_matrix <- getMutMatrix("T")
   #T_matrix <- as.matrix(T_matrix)
   gc()
   result <- rbind(A_matrix, C_matrix, G_matrix, T_matrix)
 #  result <- as.matrix(result)
   return(result)
 }# 