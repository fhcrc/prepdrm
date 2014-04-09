#' Aggregate counts per sample_name visit combination
#' @param counts Data frame, as from \code{load_counts_and_meta}
#' @param grouping Grouping values
#' @return Data frame
#' @export
aggregate_counts_per_sample <- function(counts,
                                        grouping=.(location, sample_name, visit, matched_control, reference_id, is_control, is_wt_control))
{
    agg <- function(piece) {
      data.frame(wt=sum(piece$wt),
                 mut=sum(piece$mut),
                 other=sum(piece$other),
                 non_coding=sum(piece$non_coding),
                 nonmutant=sum(piece$nonmutant),
                 n_with_mut=sum(piece$mut > 0),
                 prop_with_mut=mean(piece$mut > 0),
                 mut_counts=paste(piece$mut, collapse=','),
                 sample_ids=paste(piece$sample_id, collapse=','))
    }
    result <- ddply(counts, grouping, agg)

    ## Add proportions
    p <- result[,c('mut', 'nonmutant')]
    colnames(p) <- c('mut_prop', 'nonmutant_prop')
    p <- sweep(p, 1, rowSums(p),  '/')
    cbind(result, p)
}

#' Aggregate amino acid counts per sample
#' @param counts output of load_aa_and_meta
#' @return Data frame with counts, wide-like
#' @export
aggregate_aa_per_sample <- function(counts) {
  # Subset to direction of interest
  counts <- counts[(counts$any_indels == 'no' | counts$whole_codon_indels == 'yes') &
                   counts$direction == get_direction(counts$location), ]
  cst <- dcast(counts, location+sample_name+visit~amino_acid,
               sum,
               value.var='count',
               fill=0)
  colnames(cst)[colnames(cst) == '!'] <- 'nc'
  return(cst)
}
