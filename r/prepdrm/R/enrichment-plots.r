
#' Convert a long data frame of p-values and results to a wide one, write to a file
#' @param df Data frame. Must have columns p_value,
#' @param output_path Output path for CSV file
#' @param p_value_var Name of the p value column
#' @param cutoff Cutoff for \code{p_value_var} to consider significant
#' @param extra_vars Extra columns to include in teh output
#' @return (Invisibly) the plot
#' @export
write_wide_results <- function(df, output_path, p_value_var='p_value',
                               cutoff=0.05, extra_vars=c()) {
  cutoff_var <- sprintf('significant_%f', cutoff)
  df[, cutoff_var] <- df[, p_value_var] < cutoff
  wide_output <- reshape(df[, c('location', 'sample_name', 'visit',
                                p_value_var, extra_vars, cutoff_var,
                                'mut_prop', 'mut', 'nonmutant')],
                         direction='wide',
                         idvar=c('sample_name', 'visit'),
                         timevar='location',
                         sep='_')
  write.csv(wide_output, output_path, row.names=FALSE, na='')
  invisible(wide_output)
}

#' Visualization of mutation rates at different visits
#' @param df Data frame. Must have columns p_value,
#' @param output_path Output path for PDF
#' @param locs Locations to include
#' @param p_value_var Name of the p value column
#' @return (Invisibly) the plot
#' @export
write_enrichment_plot <- function(df, output_path,
                                  locs = c('K65R', 'M184VI', 'K70E', 'Y181CIV'),
                                  p_value_var = 'p_value') {

  df$p_value_cat <- cut(df[, p_value_var], c(0, 1e-3, 0.01, 0.05, 0.1, 1.0),
                        labels=c('<0.001', '<0.01', '<0.05', '<0.1', '>=0.1'),
                        include.lowest=TRUE)

  keep <- (!df[, 'is_control'] &
           df[, 'visit'] %in% c('SC', 'SC+1') &
           df[, 'location'] %in% locs)
  p_toplot <- df[keep, ]

  p <- ggplot(p_toplot, aes(x=location, y=sample_name, fill=p_value_cat)) +
    geom_tile(aes(color=p_value_cat)) +
    geom_text(aes(label=sprintf('%.2f%%', mut_prop * 100))) +
    facet_wrap(~visit) +
    scale_fill_manual("p-value", values=rev(brewer.pal(5,"PuBu"))) +
    scale_color_manual("p-value", values=rev(brewer.pal(5,"PuBu"))) +
    ylab('Subject') +
    xlab('Location') +
    ggtitle("Mutation percentages") +
    theme_bw(16)

  n_samples <- length(unique(p_toplot$sample_name))

  ggsave(output_path, p, width=length(locs) * 12 / 4, height=n_samples * 12/50)
  invisible(p)
}
