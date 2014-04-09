#' Adjusts p-values
#'
#' \code{adjust_p_values} adjusts p-values, by default using FDR correction,
#' for a data frame testing for enriched sites from PrEP.
#'
#' By default, only control p-values are shown adjusted, though the correction uses all.
#' Primary sites (defined in \code{prep_primary_sites}) are corrected using p-values only from those sites.
#' Secondary sites are corrected using all p-values.
#'
#' @param df Data frame
#' @param primary_separate Correct primary sites separately?
#' @param controls_only Only display corrected p-values for controls?
#' @param input_column Column containing uncorrected p-value
#' @param method passed to \code{p.adjust}
#' @return Corrected p-value
#' @export
adjust_p_values <- function(df,
                            primary_separate = TRUE,
                            controls_only = TRUE,
                            input_column = 'p_value',
                            method = 'fdr') {
  for(column in c(input_column, 'is_wt_control', 'location', 'is_control'))
    assert_that(column %in% colnames(df))

  p <- df[[input_column]]
  p_adj <- rep(NA, length(p))
  #sel <- !df[['is_wt_control']]
  sel <- seq_along(p_adj)

  p_adj[sel] <- p.adjust(p[sel], method='fdr')

  if (primary_separate) {
    sel <- sel & df[['location']] %in% prep_primary_sites
    p_adj[sel] <- p.adjust(p[sel], method='fdr')
  }

  if (controls_only) {
    # Wipe-out the adjusted p-values for non-controls
    is_control <- df[['is_control']]
    p_adj[!is_control] <- NA
  }

  p_adj
}
