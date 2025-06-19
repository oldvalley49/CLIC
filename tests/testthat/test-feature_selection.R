library(testthat)
library(Matrix)
library(devtools)

# load CLIC
devtools::load_all()



# ---------- OUTPUT CHECKS ----------
test_that("Returns correct features with unnormalized Seurat object",{
    out <- FindCLICFeatures(rna, species="human", verbose=FALSE)
    expect_setequal(out$use_features, use_features_truth)
})
rna <- NormalizeData(rna)
test_that("Returns correct features with normalized Seurat object",{
    out <- FindCLICFeatures(rna, species="human", verbose=FALSE)
    expect_setequal(out$use_features, use_features_truth)
})
test_that("Returns correct features with raw count matrix",{
    out <- FindCLICFeatures(rna_counts, species="human", verbose=FALSE)
    expect_setequal(out$use_features, use_features_truth)
})
test_that("Returns correct number of features from matrix input", {
  out <- FindCLICFeatures(mat, species="human", verbose=FALSE, nfeatures = 500)
  expect_equal(length(out$use_features), 500)
})
test_that("Returns correct number of features from Seurat object", {
  seurat <- CreateSeuratObject(counts = as(mat, 'dgCMatrix'))
  out <- FindCLICFeatures(seurat, species="human", verbose=FALSE, nfeatures = 500)
  expect_equal(length(out$use_features), 500)
})

# ---------- ERROR CHECKS ----------

test_that("Fails if nfeatures > initial_variable_features_num", {
  expect_error(
    FindCLICFeatures(mat, species="human", verbose=FALSE, nfeatures = 2000, initial_variable_features_num = 1000),
    "`nfeatures` must be less than or equal to `initial_variable_features_num`"
  )
})

test_that("Fails if initial_variable_features_num > number of rows", {
  expect_error(
    FindCLICFeatures(mat, species="human", verbose=FALSE, nfeatures = 2000, initial_variable_features_num = 20000),
    "cannot exceed the number of features"
  )
})

test_that("Fails if matrix has fewer rows than nfeatures", {
  expect_error(
    FindCLICFeatures(mat[1:100, ], species="human", verbose=FALSE, nfeatures = 20000),
    "nfeatures should be smaller than the number of features"
  )
})

test_that("Fails if species file not found", {
  expect_error(
    FindCLICFeatures(mat, species = "alien"),
    "CLIC scores for the species are not currently available"
  )
})


# ---------- WARNING CHECKS ----------

test_that("Warns missing feature names", {
  mat2 <- mat
  rownames(mat2) <- NULL
  expect_warning(
    FindCLICFeatures(mat2, species="human", verbose=FALSE),
    "No feature names were provided"
  )
})

test_that("Warns on duplicate feature names", {
  mat3 <- mat
  rownames(mat3)[1:2] <- "dup-gene"
  expect_warning(
    FindCLICFeatures(mat3, species="human", verbose=FALSE),
    "duplicate feature names"
  )
})

test_that("Warns on duplicate cell names", {
  mat4 <- mat
  colnames(mat4)[1:2] <- "dup-cell"
  expect_warning(
    FindCLICFeatures(mat4, species="human", verbose=FALSE),
    "duplicate cell names"
  )
})

