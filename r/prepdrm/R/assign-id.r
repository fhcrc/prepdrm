#' Assign a unique subject ID using first plate, barcode for each subject
#'
#' @param sample_names Sample names - for grouping
#' @param sample_ids Sample IDs, i.e. p1d1bc241
#' @return a new sample_ids
#' @export
anonymize_ids <- function(sample_ids, sample_names) {
  assert_that(is.character(sample_names) || is.factor(sample_names))
  assert_that(is.character(sample_ids) || is.factor(sample_ids))
  assert_that(length(sample_ids) == length(sample_names))

  sample_ids <- as.character(sample_ids)
  pat <- 'p(\\d)d.bc(\\d{3}).*'

  assert_that(all(grepl(pat, sample_ids)))

  sample_names <- factor(sample_names)
  lu <- tapply(sample_ids, sample_names,
               function(x) sub(pat, '\\1-\\2', sort(x)[1]))
  result <- as.vector(lu[sample_names])
  unname(result)
}
