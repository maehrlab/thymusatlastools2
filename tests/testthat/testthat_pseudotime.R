testthat::context("pseudotime")
expect_runs = function( expr ){
  expect_false( identical("fail", tryCatch( expr, finally = function() "fail") ) )
}
library(magrittr)

# load a minimal example data set (subset of thymus atlas)
thymus_test = readRDS("../testdata/thymus_test.Rdata")
thymus_test %<>% (Seurat::UpdateSeuratObject)

# For interactive use:
# thymus_test = readRDS("tests/testdata/thymus_test.Rdata")
# thymus_test %<>% (Seurat::UpdateSeuratObject)


testthat::test_that( "ProjectCells runs and gets close to recovering actual values", {
  thymus_test %<>% RunPCA(pcs.compute = 25)
  thymus_test %<>% ProjectCells(thymus_test, to_project = "tSNE_1", regressout = "nUMI")
  expect_gt(cor(thymus_test@meta.data$tSNE_1, thymus_test@dr$tsne@cell.embeddings[, 1]), 0.9)
})

testthat::test_that( "KNN projection works", {
  thymus_test %<>% master_pt(results_path = "../testdata/tmp")
  assertthat::has_name( thymus_test@meta.data, "pseudotime" )
})


testthat::test_that( "smoothing  and clustering works", {
  current_rp = "../testdata/tmp"
  thymus_test %<>% master_pt(results_path = current_rp, method = "DPT")
  gene_stats = GetDynamicGenes( thymus_test, num_periods_initial_screen = 10, pt.use = "branch_viz_1" )
  smoothers = SmoothGenes( thymus_test, genes.use = gene_stats$gene %>% head(20), pt.use = "branch_viz_1" )
  cluster_mod_etc = ClusterSmoothedGenes( thymus_test, results_path = current_rp, pt.use = "branch_viz_1",
                                          smoothers = smoothers, abcissae_kmeans = 20 )
  FacetPlotGeneClusters(
    thymus_test, results_path = current_rp, pt.use = "branch_viz_1",
    cluster_mod = cluster_mod_etc$cluster_mod,
    kmeans_features = cluster_mod_etc$kmeans_features,
    abcissae_kmeans = cluster_mod_etc$abcissae_kmeans
  )
  HeatmapGeneClusters(
    thymus_test, results_path = current_rp, pt.use = "branch_viz_1", 
    cluster_mod = cluster_mod_etc$cluster_mod,
    smoothers = smoothers, 
    genes_use = names(smoothers)
  )
})



unlink("../testdata/tmp")
