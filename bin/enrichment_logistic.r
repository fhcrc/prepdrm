#!/usr/bin/env Rscript

library(assertthat)
library(ggplot2)
library(glmnet)
library(gridExtra)
library(plyr)

## Dev package
library(methods)
library(devtools)
suppressPackageStartupMessages(load_all('r/prepdrm', export_all=FALSE))

## Functions
sample_visit <- function(x, y) {
    assert_that(!is.null(x))
    assert_that(!is.null(y))
    paste(x, y, sep='-')
}

test_sample_enrichment <- function(piece) {
    ## Response (proportion mutated) is a function of the sample visit
    ## and predicted bias from plate / genotype only.
    ## We're thus interested in the significance of the intercept term
  m0 <- glm(prop_mutated~offset(pred_log_odds) + 0,
            data=piece,
            weights=total,
            family=binomial(link='logit'))
  m1 <- glm(prop_mutated~offset(pred_log_odds),
            data=piece,
            weights=total,
            family=binomial(link='logit'))
  a <- anova(m0, m1, test='LRT')

  result <- as.data.frame(summary(m1)$coef[, c('Estimate', 'Pr(>|z|)'), drop=FALSE])
  colnames(result) <- c('log_odds', 'logit_p_value')
  result$converged <- m1$converged
  result$lrt_p_value <- ifelse(result$log_odds > 0, a[2, 'Pr(>Chi)'], 1)
  result$offsets <- paste(piece$pred_log_odds, collapse=' ')
  result
}

## Test for enrichment for a single location at a time.
enrichment_per_location <- function(df, pseudocount=0) {
    ## We'll decorate the per-sample mutant counts with assessments of significance
    sample_counts <- aggregate_counts_per_sample(df)
    sample_counts$sample_visit <- sample_visit(sample_counts$sample_name,
                                               sample_counts$visit)


    ## Fit a model from the wild-type controls
    wt_formula <- formula(prop_mutated ~ plate)
    if(grepl('^K(65|70)', df$location[1]))
        wt_formula <- update(wt_formula, ~ . + matched_control)

    df_orig <- df

    ## Add a single mutant pseudocount
    control_df <- transform(df[df$is_wt_control, ],
                            prop_mutated=(mut + pseudocount) / (total + pseudocount))

    if(sum(control_df$mut > 0) <= 1)
      return(NULL)

    X <- model.matrix(wt_formula,
                      data=control_df,
                      na.action=NULL)
    y <- cbind(1 - control_df[['prop_mutated']], control_df[['prop_mutated']])
    w <- control_df[['total']]
    #save(df, control_df, wt_formula, sample_counts, file='model.Rdata')
    glmnet.control(pmin=1e-14)
    wt_model <- cv.glmnet(x=X, y=y,
                          weights=w,
                          family='binomial')
    wt_coef <- coef(wt_model)[, 1]
    nonzero <- names(wt_coef)[wt_coef != 0]

    ## Calculate a predicted log-odds for each row
    df$pred_log_odds <- predict(wt_model, newx=model.matrix(wt_formula, df, na.action=NULL))[, 1]

    p <- ggplot(df, aes(x=paste(sample_name, visit), y=pred_log_odds)) +
      geom_point(size=1.2) +
      ggtitle(sprintf('%s\nincluded terms: %s',
                      df$location[1], paste(nonzero, collapse=', '))) +
      facet_wrap(~plate, scales='free_x', nrow=1) +
      theme(axis.text.x=element_text(angle=90, size=rel(0.7)),
            legend.position='bottom') +
      xlab('Sample / Visit') +
      ylab('Predicted log-odds')
    if(grepl('^K(65|70)', df$location[1]))
      p <- p + geom_point(aes(color=matched_control), size=1.2)

    ## Fit a model per subject-sample
    estimates <- ddply(df,
                       .(sample_visit),
                       test_sample_enrichment)

    ## One-sided test - ignore any samples with log_odds < 0 (*less likely*
    ## than controls to have mutation), divide p-value by 2
    estimates[['lrt_p_value']] <- ifelse(estimates[['log_odds']] <= 0,
                                         1,
                                         estimates[['lrt_p_value']] / 2)
    estimates$lrt_signif_0.05 <- estimates[['lrt_p_value']] < 0.05 & estimates[['log_odds']] > 0
    estimates$odds_ratio <- exp(estimates[['log_odds']])

    local({
      m <- merge(df[, c('sample_visit', 'pred_log_odds')],
                  estimates[, c('sample_visit', 'lrt_signif_0.05', 'log_odds')],
                  all=TRUE)
      p2 <- ggplot(m, aes(x=pred_log_odds, y=log_odds)) +
        theme(legend.position='bottom') +
        ggtitle('Predicted vs. observed') +
        geom_hline(yintercept=0, linetype='dashed') +
        xlab('Predicted log-odds') +
        ylab('Fit log-odds')
      if(length(unique(m$pred_log_odds)) == 1)
        p2 <- p2 + geom_boxplot(aes(x=factor(pred_log_odds), fill=lrt_signif_0.05))
      else
        p2 <- p2 + geom_point(aes(color=lrt_signif_0.05)) +
          geom_abline(slope=1)

      grid.arrange(p, p2, ncol=2, widths=c(9, 2))
    })

    m <- merge(sample_counts, estimates, by='sample_visit', all.x=TRUE)
    m[, which(colnames(m) != 'sample_visit')]
}


## Script stuff
main <- function(sample_meta_path, counts_path, outfile, outpdf, pseudocount=0) {
    counts <- load_counts_and_meta(counts_path, sample_meta_path)
    #counts <- counts[grepl('K103|K65|M184|Y181', counts[['location']]), ]

    counts$sample_visit <- sample_visit(counts$sample_name, counts$visit)
    factor_cols <- c('plate', 'matched_control', 'sample_visit', 'sample_id', 'visit', 'sample_name')
    counts$plate <- factor(counts$plate)
    counts[, factor_cols] <- lapply(counts[, factor_cols], as.factor)
    counts[, 'matched_control'] <- relevel(counts[, 'matched_control'], 'MB2059_pol')

    counts$total <- counts$mut + counts$nonmutant
    counts$prop_mutated <- counts$mut / counts$total

    #counts <- subset(counts, location %in% c('K103NS', 'Y181C', 'M184VI', 'K65R', 'K65N'))

    message('testing enrichment at ', length(unique(counts$location)), ' locations.')
    pdf('logit_model.pdf', width=30, useDingbats=FALSE)
    theme_set(theme_bw(12))
    results <- ddply(counts, .(location), enrichment_per_location,
                     .progress='text',
                     pseudocount=pseudocount)
    dev.off()

    # FDR correct for *controls only*
    results$lrt_p_adj <- adjust_p_values(results, input_column = 'lrt_p_value',
                                         controls_only = FALSE)

    write.csv(results, outfile, row.names=FALSE, na='')
    write_enrichment_plot(results, outpdf, p_value_var='lrt_p_value')
    invisible(results)
}


if(!interactive()) {
  args <- commandArgs(TRUE)
  stopifnot(length(args) %in% 4:5)
  sample_meta_path <- args[1]
  counts_path <- args[2]
  outfile <- args[3]
  outpdf <- args[4]

  pseudocount <- 0
  if(length(args) == 5) pseudocount <- as.integer(args[5])

  # Run it
  main(sample_meta_path, counts_path, outfile, outpdf, pseudocount=pseudocount)
}
