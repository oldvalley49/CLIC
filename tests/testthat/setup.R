# load testing data (BMMC-s4d8) for ground truth
load(testthat::test_path("testdata/BMMC-s4d8.Rdata"))
rna_counts <- GetAssayData(rna, layer='counts')
use_features_truth <- read.csv(testthat::test_path("testdata/BMMC-s4d8_clic_5000_features_truth.csv"))
use_features_truth <- use_features_truth$gene

# Create public simulated data
set.seed(42)
mat <- matrix(rpois(10000000, lambda = 10), nrow = 10000)
rownames(mat) <- paste0("gene", 1:10000)
colnames(mat) <- paste0("cell", 1:1000)