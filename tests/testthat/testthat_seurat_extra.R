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


testthat::test_that( "AvailableData runs and gives back genes and essentials in a bare-bones object", {
  should_have = c( "nGene"    ,     "nUMI"      ,    "orig.ident" , "ident", "eday",
                   rownames(thymus_test@data))
  expect_true(all(should_have %in% AvailableData(thymus_test)))
})


testthat::test_that( "FetchDataZeroPad works for things we have and things we lack", {
  have = c("nGene"    ,     "Foxn1" )
  lack = "iq34gobafrovurbryeubv78b3"
  expect_warning(       FetchDataZeroPad( thymus_test,          lack ) )
  X = suppressWarnings( FetchDataZeroPad( thymus_test, c( have, lack ) ) )
  expect_equal( X[1:2], FetchData( thymus_test, have ) )
  expect_equal( X[[3]]  , rep( 0, nrow( X ) ) )
  expect_equal( 3, ncol( X ) )
})


testthat::test_that( "FetchDataZeroPad works for things we have and things we lack", {
  queries = c( "nGene"    ,     "Foxn1"      ,    "iq34gobafrovurbryeubv78b3" )
  expect_warning(      FetchDataZeroPad(thymus_test, queries))
  X = suppressWarnings(FetchDataZeroPad(thymus_test, queries))
  atat(all(0==X$iq34gobafrovurbryeubv78b3))
  atae(ncol(X), 3)
})


testthat::test_that( "SubsetDataFlex works", {
  thymus_test %<>% SubsetDataFlex(vars.use = "eday", predicate = "eday==12.5")
  expect_true(thymus_test@meta.data$eday %>% table %>% names %>% equals(12.5))
})

testthat::test_that( "FindMarkersFlex can match FindMarkers", {
  gu = thymus_test@data %>% rownames %>% sample(100)
  suppressWarnings({ X = FindMarkersFlex(object = thymus_test,
                                         ident.use = "ident", test.use = "wilcox",
                                         ident.1 = "C1", ident.2 = "C2", genes.use = gu ) })
  suppressWarnings({ Y = FindMarkers    (thymus_test,
                      ident.1 = "C1", ident.2 = "C2", genes.use = gu ) })
  Y = Y[order(-Y$avg_logFC), ]
  expect_equal(X[1:5], Y[1:5])
})


testthat::test_that( "SeuratPie runs", {
  expect_runs(SeuratPie(thymus_test, ident.use = "ident" ))
})



testthat::test_that( "get_similar_genes runs ", {
  expect_runs({
    get_similar_genes(thymus_test, "Epcam", 10 )
  })
})


testthat::test_that( "TACS runs with and without faceting and cutoffs", {
  expect_runs({
    TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc" )
    TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday" )
    TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday", cutoffs = c(0.7, 0.7) )
    TACS(thymus_test, gene1 = "Epcam", gene2 = "Ptprc", facet_by = "eday", cutoffs = c(0.7, 0.7), density = T )
  })
})

testthat::test_that( "", {
  best_params = expand.grid(excess_var_cutoff = c( 0.5 ),
                            log_expr_cutoff = c( 0.1 ), 
                            num_pc = c( 5 ), 
                            clust_method = "SNN",
                            clust_granularities_as_string = "0.5",
                            stringsAsFactors = F)
  # Bare bones
  thymus_test %<>% explore_embeddings( results_path = "test_results",
                                       regress_out = c(),
                                       all_params = best_params[1, ],
                                       plot_each_iteration = c("Foxn1", "orig.ident", "eday") )
  GetDimMarchenkoPastur(thymus_test)
  thymus_test %<>% LabelDoublets(rep_field = "species")
  
  unlink("./test_results", recursive=TRUE )
})

# To clean up by hand:
# unlink("~/Desktop/software_projects/thymusatlastools2/tests/test_results", recursive=TRUE )


