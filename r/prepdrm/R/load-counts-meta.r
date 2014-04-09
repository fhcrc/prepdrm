
#' Determines whether a mutation should be read by the forward or reverse read
#' @param s Mutation string
#' @return 'forward' or 'reverse'
#' @export
get_direction <- function(s) {
    s <- as.integer(sub('^[^0-9]+(\\d+).*$', '\\1', s))
    ifelse(s <= 140, 'forward', 'reverse')
}


#' Columns to extract from metadata
sample_cols <- c('plate', 'barcode_id', 'reference_id', 'is_control',
                 'is_wt_control', 'matched_control', 'direction',
                 'sample_name', 'sample_id', 'visit')

#' Default metadata path
default_meta_path <- 'output/prep_merged_metadata.csv'

#' Load and merge mutation counts and metadata
#' @param count_path Path to counts file
#' @param meta_path Path to sample metadata
#' @return Data frame with merged results
#' @export
load_counts_and_meta <- function(count_path = 'output/prep_merged_classified_counts.csv',
                                 meta_path = default_meta_path) {
    assert_that(is.string(count_path))
    assert_that(is.readable(count_path))
    counts <- read.csv(count_path, as.is=TRUE)
    assert_that(is.string(meta_path))
    assert_that(is.readable(meta_path))
    sample_meta <- read.csv(meta_path, as.is=TRUE)

    m <- merge(counts,
               sample_meta[, sample_cols],
               by.x=c('plate', 'barcode', 'direction'),
               by.y=c('plate', 'barcode_id', 'direction'),
               all.x=TRUE)
    m[, 'nonmutant'] <- m[, 'wt'] + m[, 'other']
    m <- m[m$direction == get_direction(m$location), ]
    m
}

#' Load and merge amino acid counts and metadata
#' @param count_path path to counts
#' @param meta_path path to metadata csv
#' @return Data frame with merged results
#' @export
load_aa_and_meta <- function(count_path = 'output/prep_merged_classified.csv',
                             meta_path = default_meta_path) {
    assert_that(is.string(count_path))
    assert_that(is.readable(count_path))
    counts <- read.csv(count_path, as.is=TRUE)
    assert_that(is.string(meta_path))
    assert_that(is.readable(meta_path))
    sample_meta <- read.csv(meta_path, as.is=TRUE)

    m <- merge(counts,
               sample_meta[, sample_cols],
               by.x=c('sample'),
               by.y=c('sample_id'),
               all.x=TRUE)
    return(m)
}
