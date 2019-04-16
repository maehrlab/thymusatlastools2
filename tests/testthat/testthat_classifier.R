testthat::context("seurat extras")
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


testthat::test_that( "Classifier runs and gives back prediction in metadata.", {
  expect_true(  is.null(  thymus_test@meta.data$classifier_ident ) )
  expect_runs( {
    mlr_mod = TrainClassifierMLR( thymus_test,
                                     genes.use = get_mouse_tfs(),
                                     do_knn_smooth = T, 
                                     do.save = F ) 
  } )
  expect_runs( {
    mlr_mod = TrainClassifierMLR( thymus_test,
                                     genes.use = get_mouse_tfs(),
                                     do_knn_smooth = F, 
                                     do.save = F )
  } )
  expect_runs( get_classifier_coefs(mlr_mod) )
  expect_runs( plot_classifier_coefs(mlr_mod) )
  dir.create("mlr_test_temp")
  expect_runs( DisplayTrainingResults(dge_train = thymus_test,
                                      results_path = "mlr_test_temp",
                                      mlr_mod,
                                      ident.use = "ident") )
  file.remove( "mlr_test_temp", recursive = T )
  expect_runs( { thymus_test = classify_mlr( thymus_test, mlr_mod ) } )
  expect_runs( { thymus_test = classify_mlr( thymus_test, mlr_mod, do_knn_smooth = T, k = 3 ) } )
  expect_false(  is.null(  thymus_test@meta.data$classifier_ident ) )
  expect_runs(  Seurat::FetchData( thymus_test, "classifier_probs_C" %>% paste0(1:10) ) )
  expect_runs(  all(paste0("C", 1:10) %in% Seurat::FetchData( thymus_test, "classifier_ident" ) ) )
})

