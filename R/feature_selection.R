# internal function
# var_genes and cor_genes must be ordered
intersect_until_n <- function(object, var_genes, CLIC_genes, n) {

   # initialize an empty vector to store the result
   result <- c()        
   clic_num = 0
   
   for (gene in CLIC_genes) {
      clic_num = clic_num + 1
      if (gene %in% var_genes) {
         result <- c(result, gene)
      }
      if (length(result) == n) {
         return(list(use_features = result, clic_num=clic_num))
      }
      # border case, if genes with CLIC scores available are not enough to get n features
      # this may happen if the number of computed variable genes is too close to n
      if (clic_num == length(CLIC_genes)) {
        remainder = var_genes[!(var_genes %in% CLIC_genes)]
        num_to_be_filled = n - length(result)
        result <- c(result, remainder[1:num_to_be_filled])
        return(list(use_features = result, clic_num=clic_num))
      }
   }
   
}


#' Select features for scRNA-seq + scATAC-seq integration using CLIC scores
#'
#' This function identifies high-confidence features for cross-modal integration
#' by combining traditional highly variable feature (HVF) selection with
#' CLIC scores, which quantify the confidence of feature-to-feature links.
#'
#' In the context of scRNA-seq and scATAC-seq integration, the CLIC score is computed
#' as the empirical correlation between gene expression (from RNA) and gene activity
#' (from ATAC) for each gene. Refer to our publication for details on how the CLIC score is computed.
#'
#' The input can be either a Seurat object (recommended) or a raw counts expression matrix.
#' If a matrix is provided, it will be converted into a Seurat object internally.
#'
#' @param object A Seurat object or a gene expression counts matrix.
#' @param score_name A character string specifying the type of CLIC score used depending on species and gene activity computation.
#'    Currently, "human_signac" and "mouse_signac" are supported.
#' @param nfeatures Number of features to select, equivalent to `nfeatures` in Seurat's \code{FindVariableFeatures()}. 
#'    Default is 2000.
#' @param initial_variable_features_num Number of variable features to compute before applying CLIC selection. 
#'    A higher value gives more weight to the CLIC scores. Default is 5000.
#' @param verbose Logical. Whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
#' \describe{
#'   \item{use_features}{Character vector of selected features.}
#'   \item{clic_num}{Integer indicating how many CLIC-ranked features were required to reach \code{nfeatures}. 
#'    The intersection of the top \code{initial_variable_features_num} variable genes and the top \code{clic_num} 
#'    CLIC-ranked genes determines the final feature set.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- FindCLICFeatures(object = obj, score_name = "human")
#' use_features <- out$use_features
#' }
#'
#' @author Tomoya Furutani \email{furutomo.49@@gmail.com}
#' 
#' @importFrom Seurat CreateSeuratObject NormalizeData FindVariableFeatures HVFInfo GetAssayData SetAssayData VariableFeatures
#' @importFrom methods as
#' @importClassesFrom Matrix Matrix
#' 
#' @export
FindCLICFeatures <- function(object, score_name, nfeatures=2000, initial_variable_features_num=5000, verbose=TRUE) {

   # check data properties

   if (nrow(object) < nfeatures) {
      stop("nfeatures should be smaller than the number of features in expression matrix")
   }
   if (is.null(rownames(object))) {
      rownames(object) <- c(1:nrow(object))
      warning("No feature names were provided, outputting the row numbers of selected features")
   }

   if (is.null(colnames(object))) {
      colnames(object) <- c(1:ncol(object))
   } else if (length(unique(colnames(object))) < ncol(object)) {
      warning("There are duplicate cell names, making cell names unique")
      colnames(object) <- make.unique(colnames(object))
   } else if (length(unique(rownames(object))) < nrow(object)) {
      warning("There are duplicate feature names, making feature names unique")
      rownames(object) <- make.unique(rownames(object))
   }

   # check parameters

   if (nfeatures > initial_variable_features_num) {
      warning('number of features in dataset were less than `initial_variable_features_num`')
   }
   
   scores_fp <- system.file('extdata', paste0(score_name,'.csv'), package='CLIC')
   if (!file.exists(scores_fp)){
      stop('CLIC scores with the name are not currently available')
   } else {
      # load CLIC scores here
      scores <- read.csv(scores_fp)
   }

   # check input: SeuratObject, matrix, or other (raise error)
   if (inherits(object, "Seurat")) {

      # check if RNA is an assay
      if (!"RNA" %in% names(object@assays)) {
         stop("Default RNA assay not found")
      }
      # check if already normalized, if not, normalize the data
      data_mat <- GetAssayData(object, layer='data')
      if (nrow(data_mat) == 0 || ncol(data_mat) == 0) {
         if (verbose) message("Seurat object provided is not normalized. Normalizing...")
         object <- NormalizeData(object)
      } else {
         if (verbose) message("Seurat object appears to be already normalized. Continuing without normalization...")
      }
      rm(data_mat)
   } else if (inherits(object, 'Matrix') || inherits(object, 'matrix')) {
      if (verbose) message(paste0("Treating input matrix as counts matrix"))
      if (!inherits(object, 'dgCMatrix')){
         object <- as(object, 'dgCMatrix')
      }
      object <- CreateSeuratObject(counts=object)
      object <- NormalizeData(object)
   } else {
      stop("Input must be either a Seurat object or a matrix")
   }

   # find variable features
   object <- FindVariableFeatures(object, nfeatures=initial_variable_features_num, verbose=verbose)
   var.info <- HVFInfo(object)
   var.info.sorted <- var.info[order(var.info$variance.standardized, decreasing=TRUE), ]
   var.features.sorted <- head(rownames(var.info.sorted), n=initial_variable_features_num)
   # check if var features match
   stopifnot(all(var.features.sorted %in% VariableFeatures(object)))

   # order CLIC genes
   scores.sorted <- scores[order(scores$CLIC_score, decreasing=TRUE), ]
   CLIC_genes.sorted <- scores.sorted$gene

   # get CLIC features
   selection <- intersect_until_n(object, var.features.sorted, CLIC_genes.sorted, n=nfeatures)

   return(selection)
}
