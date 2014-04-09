context("data")

test_that("Data are accessible", {
  expect_equal('^(K65|D67)', genotype_specific_sites)
})
