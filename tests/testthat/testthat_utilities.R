testthat::context("utilities")

expect_runs = function( expr ){
  expect_false( identical("fail", tryCatch( expr, finally = function() "fail") ) )
}

testthat::test_that( "aggregate_nice works for means and sums", {
  test1 = data.frame(1:6, 2, 3) 
  rownames(test1) = letters[1:6]
  expect_true( all(as.matrix( aggregate_nice( test1, letters[1:6], mean ) == test1 )) )
  expect_equal( aggregate_nice( test1, rep("A", 6), sum), 
                matrix( c(21, 12, 18), nrow = 1, dimnames = list("A", c("X1.6", "X2", "X3"))) ) 
})


testthat::test_that( "replace_with_int_rank works for integers and characters", {
  expect_equal(replace_with_int_rank(1:5), 1:5)
  expect_equal(replace_with_int_rank(c("a", "b", "a", "c", "c")), c(1, 2, 1, 3, 3))
  expect_equal(replace_with_int_rank(c("0", "2", "0", "4", "4")), c(1, 2, 1, 3, 3))
})

testthat::test_that( "replace_with_int_rank works for integers and characters", {
  expect_equal(Capitalize(c("FOXN1", "PSMB11")), c("Foxn1", "Psmb11") )
})


testthat::test_that( "strip_suffix works for removing '.pdf'", {
  expect_equal( strip_suffix("blah.pdf", ".pdf"), "blah")
  expect_equal( strip_suffix("blah",     ".pdf"), "blah")
})


testthat::test_that( "get_preimage works for simple examples", {
  expect_equal( get_preimage( map = setNames( LETTERS, letters) ),
        as.list( setNames( letters, LETTERS) )  )
  to_invert = setNames(      c("A", "B", "DUPE", "DUPE", "DUPE2", "DUPE2"),
                             nm = c("a", "b", "c",    "d", "e",    "f") )
  inverse = list(A = "a",
                 B = "b",
                 DUPE  = c( "c", "d" ) ,
                 DUPE2 = c( "e", "f" ) )
  expect_equal( get_preimage( map = to_invert ), inverse )
})

testthat::test_that( "div_by_max , div_by_sum , and standardize work for zeroes and constants and simple examples", {
  z = rep(0, 5)
  e = c(5, 0, 0, 0, 0)
  expect_equal(div_by_max(z), z)
  expect_equal(div_by_max(e), e/5)
  expect_equal(div_by_sum(z), z)
  expect_equal(div_by_sum(e), e/5)
  expect_equal(standardize(z), z)
  expect_equal(standardize(z+1), z)
  expect_equal(var(standardize(c(-1,0,1))), 1)
})


testthat::test_that( "top_n_preserve_rownames works even when name of temp col is taken", {
  # Make sure an adversarial case -- temp column name already taken -- works out ok
  expect_equal(  top_n_preserve_rownames(x = data.frame(rownames_tempcol = 10:6), 3, rownames_tempcol),
                 data.frame(rownames_tempcol = 10:8) ) 
})

library(magrittr)
thymus_test = readRDS("../testdata/thymus_test.Rdata")
thymus_test %<>% UpdateSeuratObject

testthat::test_that( "convert_species_dge runs (m2h and h2m)", {
  expect_runs({
    thymus_test %>%
      convert_species_dge( to = "human", from = "mouse") %>% 
      convert_species_dge( to = "mouse", from = "human") 
    })
})




