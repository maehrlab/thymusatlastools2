testthat::context("data handling")

# load a minimal example data set (subset of thymus atlas)
# thymus_test = readRDS("tests/testdata/thymus_test.Rdata")
thymus_test =  readRDS("../testdata/thymus_test.Rdata")
thymus_test  = seuratify(thymus_test@raw.data)
thymus_test %<>% add_maehrlab_metadata( "eday")

expect_runs = function( expr ){
  expect_false( identical("fail", tryCatch( expr, finally = function() "fail") ) )
}


testthat::test_that( "Read10X works on v3", {
  expect_silent( Read10X( "../testdata/filtered_feature_bc_matrix" ) )
})

testthat::test_that( "Write10X inverts Read10X", {
  dir.create.nice(      "temp_write10x")
  Write10X(thymus_test, "temp_write10x", include_all = T)
  tt2 = Read10X(        "temp_write10x")
  me = magrittr::equals
  expect_true(thymus_test@raw.data %>% rownames %>% me(rownames(tt2)) %>% all)
  expect_true(thymus_test@raw.data %>% colnames %>% me(colnames(tt2)) %>% all)
  expect_true(thymus_test@raw.data %>%              me(        (tt2)) %>% all)
  unlink("temp_write10x", recursive = T)
})


testthat::test_that( "dge_merge_list works", {
  tt1 = thymus_test@raw.data
  tt2 = tt1
  colnames(tt2) = paste0(colnames(tt1), "_dupe")
  expect_runs( dge_merge_list(list(tt1, tt2)) )
})


testthat::test_that( "seuratify works", {
  tt1 = thymus_test@raw.data
  expect_runs(seuratify(tt1, results_path = "test_results"))
})
unlink("./test_results", recursive=TRUE )


expect_warning({thymusatlastools2::load_maehrlab_data("hTEP_Multiplex_H1LucGFP")})

