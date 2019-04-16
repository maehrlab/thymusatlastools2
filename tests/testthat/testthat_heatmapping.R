testthat::context("heatmapping and plotting")
expect_runs = function( expr ){
  expect_false( identical("fail", tryCatch( expr, finally = function() "fail") ) )
}
library(magrittr)

# load a minimal example data set (subset of thymus atlas)
thymus_test = readRDS("../testdata/thymus_test.Rdata")
thymus_test %<>% (Seurat::UpdateSeuratObject)


testthat::test_that( "custom_feature_plot runs", {
  expect_runs({
    custom_feature_plot(thymus_test)
    custom_feature_plot(thymus_test, "Epcam")
    custom_feature_plot(thymus_test, "Epcam", mode = "regular")
    custom_feature_plot(thymus_test, "Epcam", mode = "points")
    custom_feature_plot(thymus_test, "Epcam", mode = "overplot_adjust")
    custom_feature_plot(thymus_test, "Epcam", mode = "umidensity")
    custom_feature_plot(thymus_test, "Epcam", mode = "umi_density")
    custom_feature_plot(thymus_test, "Epcam", mode = "umi_density", n = 100, h = 0.01)
    custom_feature_plot(thymus_test, mode = "umi_density")
  })
})


testthat::test_that( "DoHeatmapFast runs", {
  # expect_runs({
  #   DoHeatmapFast(thymus_test, results_path = ".",
  #                 genes.preview = c("Actb", "Foxn1", "Krt5"),
  #                 genes.label = c("Actb", "Foxn1", "Krt5"),
  #                 ident.use = "eday")
  # })
  # file.remove("./heatmap_cellwise_PREVIEW.pdf")
  # file.remove("./heatmap_cellwise.pdf")
  expect_runs({
    make_heatmap_for_table( dge = thymus_test,
                            genes_in_order = c("Actb", "Foxn1", "Krt5"),
                            genes_to_label = c("Actb", "Foxn1", "Krt5"),
                            desired_cluster_order = (12:16) + 0.5, 
                            aggregator = mean,
                            normalize = "row",
                            norm_fun = div_by_max,  
                            ident.use = "eday" )
  })
})


