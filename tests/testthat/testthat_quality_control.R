testthat::context("quality_control_qc")


expect_runs = function( expr ){
  expect_false( identical("fail", tryCatch( expr, finally = function() "fail") ) )
}

thymus_test = readRDS("../testdata/thymus_test.Rdata")
thymus_test %<>% (Seurat::UpdateSeuratObject)

testthat::test_that("check_xist_pure works on at least one example", {
  expect_runs(check_xist_pure(raw_dge = thymus_test@raw.data, 
                              rep_name = "e12_5wholeThy_venus_rep2",
                              results_path = "."))
  unlink("e12_5wholeThy_venus_rep2_Xist_vs_Y_genes.txt")
})



