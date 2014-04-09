context("assign-id")

test_that("IDs are assigned as expected", {
  sample_ids <- c('p1d1bc241', 'p1d1bc242', 'p2d1bc242', 'p2d2bc280', 'p2d1bc280', 'p2d1bc281')
  sample_names <- c('A', 'A', 'A', 'B', 'B', 'B')
  expected <- c(rep('1-241', 3), rep('2-280', 3))
  lu <- anonymize_ids(sample_ids, sample_names)
  expect_equal(expected, lu)

  # Randomize
  set.seed(1)
  o <- sample(seq_along(sample_ids), length(sample_ids))
  expect_equal(expected[o], anonymize_ids(sample_ids[o], sample_names[o]))
})
test_that("Invalid inputs are rejected", {
  expect_error(anonymize_ids(c('A', 'A', 'A'),
                             c('p1d1bc241', 'p1d1bc242', 'p2d1bc242', 'p2d2bc280')))
})
