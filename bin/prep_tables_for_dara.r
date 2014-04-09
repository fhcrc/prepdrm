#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(RColorBrewer)
library(reshape2)
library(scales)
library(devtools)
load_all('r/prepdrm', export_all=FALSE)

# Global Style
output_width <- 8.5
theme_set(theme_bw(16) + theme(panel.grid = element_blank(),
                               strip.background = element_blank()))

# Plotting functions
sample_plot <- function(df, include_text = TRUE) {
  nd <- df[df$mut_prop_label == 'no data', ]
  done <- df[df$mut_prop_label != 'no data', ]
  p <- ggplot(df, aes(x=visit, y=subj_id, fill=mut_prop_cat)) +
    geom_tile() +
    facet_grid(treatment_arm~location, scales="free_y", space="free") +
    text_scale +
    fill_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Subject') +
    xlab('Visit')
  if (nrow(nd) > 0)
    p <- p +
      geom_text(aes(label=mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, mut_prop >= 0.1)),
                data = nd)
  if (include_text)
    return(p +
           geom_text(aes(label=mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, mut_prop >= 0.1)), data = done))
  else
    return(p)
}

sample_pvalue_plot <- function(df, include_text = TRUE) {
  df$p_value_cat <- categorize_p_values(df[, 'p_adj'])
  df$p_value_cat[is.na(df$mut_prop_cat)] <- NA
  nd <- df[df$mut_prop_label == 'no data', ]
  done <- df[df$mut_prop_label != 'no data', ]
  p <- ggplot(df, aes(x=visit, y=subj_id, fill=p_value_cat)) +
    geom_tile() +
    facet_grid(treatment_arm~location, scales="free_y", space="free") +
    p_value_fill +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Subject') +
    xlab('Visit')
  if (nrow(nd) > 0)
    p <- p +
      geom_text(aes(label = mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, mut_prop >= 0.1)),
                data = nd)
  if (include_text)
    return(p +
           geom_text(aes(label=mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, mut_prop >= 0.1)), data = done))
  else
    return(p)
}

mixture_theme <- theme(panel.grid = element_blank(),
                       strip.background = element_blank(),
                       strip.text = element_blank())
mixture_plot <- function(df) {
  ggplot(df, aes(x = location, y = mixture_perc_expected, fill = mut_prop_cat)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100))) +
    fill_scale +
    mixture_theme +
    ylab('Mixture') +
    xlab('Location') +
    ggtitle(paste("Genotype", df$geno_num[1]))
}

mixture_pvalue_plot <- function(df) {
  df$p_value_cat <- categorize_p_values(df[, 'p_adj'])
  ggplot(df, aes(x = location, y = mixture_perc_expected, fill = p_value_cat)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100), color = p_adj < 0.01)) +
    p_value_fill +
    text_scale +
    mixture_theme +
    ylab('Mixture') +
    xlab('Location') +
    ggtitle(paste("Genotype", df$geno_num[1]))
}

mixture_level_plot <- function(df) {
  cst <- dcast(df, sample_name ~ location, value.var='mut_prop')
  lm <- lm(M184VI ~ K65R + 0, cst)
  ct <- cor.test(cst$K65R, cst$M184VI)
  ggplot(cst, aes(x = K65R, y = M184VI)) +
    geom_abline(slope=1, linetype='dashed', color='#444444') +
    geom_point() +
    coord_equal() +
    scale_x_continuous(labels=percent, limits=c(0, 0.22)) +
    scale_y_continuous(labels=percent, limits=c(0, 0.22))
}

mixture_singlepos_plot <- function(df) {
  ggplot(df, aes(x = location, y = mixture_perc_expected, fill = p_adj < 0.05 & n_with_mut == 1)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100), color = p_adj < 0.05 & n_with_mut == 1)) +
    scale_fill_manual(values=c('white', brewer.pal(3, 'Blues')[3]), guide=FALSE) +
    text_scale +
    mixture_theme +
    ylab('Mixture') +
    xlab('Location') +
    ggtitle(paste("Genotype", df$geno_num[1]))
}

wt_plot <- function(df) {
  ggplot(df, aes(x = location, y = subj_id, fill = p_adj < 0.05)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100), color = p_adj < 0.05)) +
    scale_fill_manual(values=c(bd_color, brewer.pal(3, 'Blues')[3]), guide = FALSE) +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Sample') +
    xlab('Location')
}

wt_table <- function(df, f = function(x) sprintf('%.2f%%', mean(x * 100))) {
  sample <- gsub('(^\\d-|_?wt|_pol|-\\d$|MB)', '', df$sample_name)
  df$genotype <- factor(sample,
                        levels = paste('2059', c('', '_F', '_E', '_G', '_H'), sep=''),
                        labels = c('2059', 'F', 'E', 'G', 'H'))
  dcast(melt(df[, c('location', 'genotype', 'mut_prop')], measure.vars=3),
        location~genotype, f)
}

wt_error_by_site <- function(df) {

  df$genotype <- factor(df$group,
                        levels = paste('MB2059_pol', c('', '_G', '_H', '_E', '_F'), sep=''),
                        labels = c('1', '2', '3', '4', '5'))
    #geom_point(position='jitter') +
  p <- ggplot(df, aes(x=genotype, y=mut_prop, fill=genotype)) +
    geom_boxplot() +
    facet_wrap(~location, nrow=1) +
    theme(legend.position='none') +
    scale_y_continuous(labels=percent) +
    xlab("Genotype at K65") +
    ylab("Reads classified as mutant")
}


########################
## Load and plot data ##
########################

load_data <- function() {
  p_values <- transform(read.csv('output/prep_merged_fisher.csv', as.is=TRUE),
                        subj_id = anonymize_ids(sample_ids, sample_name))
  p_values <- transform(p_values, subj_id = factor(subj_id, levels=rev(sort(levels(subj_id)))))

  # Treatment arm, expected mut.
  treatment_arm <- read.csv('metadata/treatment_arm.csv')
  p_values$treatment_arm <- factor(treatment_arm[match(p_values[['sample_name']], treatment_arm[['sample_name']]), 'arm'],
                                   levels = c('FTC/TDF', 'TDF', 'Placebo'))
  p_values <- join(p_values, read.csv('metadata/mixture_proportions.csv', as.is = TRUE),
                   by='sample_name', match = 'first')
  p_values$mixture_perc_expected <- ordered(p_values$perc_expect,
                                            levels = sort(unique(p_values$perc_expect)),
                                            labels = paste(sort(unique(p_values$perc_expect)), '%', sep=''))

  # Add some annotations
  p_values <- mutate(p_values,
                     # Apply cutoff
                     above_detection = p_adj <= 0.05,
                     # Bin the proportion mutated, with no level for those above detection
                     mut_prop_cat = cut(mut_prop,
                                        c(0, .01, .05, .10, .20, 1.00),
                                        labels=c('<= 1%', '1 - 5%', '5 - 10%', '10 - 20%', '> 20%'),
                                        include.lowest=TRUE),
                     # A % label for above-detection items
                     mut_prop_label = ifelse(above_detection, sprintf('%.2f%%', mut_prop * 100), ''),
                     # Recode SC+3 as SC+1, recode SC+1repeat as SC+1
                     visit = sub(' ?repeat', '', ifelse(visit == 'SC+3', 'SC+1', visit)),
                     # Order locations by position
                     location = reorder(location, as.integer(sub('^[A-Z]+(\\d+).*$', '\\1', as.character(location))), max))

  # Add bd category
  p_values <- mutate(p_values,
                     mut_prop_cat = ordered(ifelse(!above_detection, 'bd',
                                                   as.character(mut_prop_cat)),
                                            levels=c('bd', levels(mut_prop_cat))))

  # Order samples by maximum mutation percent across primary sites
  prim <- p_values[p_values[['location']] %in% prep_primary_sites, ]
  l <- levels(reorder(prim[['subj_id']], ifelse(prim[['above_detection']], prim[['mut_prop']], 0), max))
  p_values$subj_id <- factor(p_values$subj_id, levels=l)

  return(p_values)
}

annotate_missing_visits <- function(df) {
  ddply(df, .(location, subj_id), function(piece) {
    for(v in c('SC', 'SC+1')) {
      if(!(v %in% piece$visit)) {
        # Add row
        piece <- rbind(piece, mutate(piece[1, ], visit=v, mut_prop_label='no data', mut_prop_cat=NA, above_detection=NA))
      }
    }
    return(piece)
  })
}

p_values <- load_data()

# Subset to samples, sites of interest
to_plot <- subset(p_values,
                  !is_control & visit %in% c('SC', 'SC+1') &
                  location %in% prep_primary_sites)
# Add rows with missing values for missing samples
to_plot <- annotate_missing_visits(to_plot)

# Equivalent for mixture controls
wt_to_plot <- subset(p_values, is_wt_control & location %in% prep_primary_sites)

# Number of subjects
n_subjects <- length(unique(to_plot$subj_id))

######################
# Big by subject table
######################
pdf('3_prep_table_all.pdf', width=output_width, n_subjects * 12/50)
tryCatch(print(sample_plot(to_plot, include_text = TRUE)), finally=dev.off())
write.csv(to_plot, '3_prep_table_all.csv', row.names=FALSE)

n_splits <- 2
pdf('prep_table_all_split.pdf', width=output_width, n_subjects * 12/50 / n_splits)
tryCatch({
  sp <- interaction(to_plot$subj_id, to_plot$treatment_arm,  drop=TRUE)
  pv <- mutate(to_plot, split=cut(as.integer(sp), breaks=n_splits))
  d_ply(pv, .(split), function(piece) {
    print(sample_plot(piece))
  })
}, finally=dev.off())

d_ply(to_plot, .(treatment_arm), function(piece) {
  arm <- sub('/', '-', piece$treatment_arm[1])
  n_subjects <- length(unique(piece$subj_id))
  write.csv(piece, sprintf('5_prep_table_all_%s.csv', arm),
            row.names=FALSE)
  pdf(sprintf('4_prep_table_all_%s.pdf', arm),
      width = output_width,
      height = n_subjects * 12/50)
  tryCatch({
    print(sample_plot(piece))
  }, finally=dev.off())
})


#################
# Mixture plots
#################
local({
  p <- function(d, prefix, width=output_width * 0.66) {
    h <- output_width * length(unique(d$subj_id)) / 50
    pdf(paste(prefix, '_p.pdf', sep=''), width=width, height=h)
    tryCatch({
      d_ply(d, .(geno_num), function(p) print(mixture_pvalue_plot(p)))
    }, finally=dev.off())
    pdf(paste(prefix, '.pdf', sep=''), width=width, height=h)
    tryCatch({
      d_ply(d, .(geno_num), function(p) print(mixture_plot(p)))
      d_ply(d, .(geno_num), function(p) print(mixture_singlepos_plot(p)))
    }, finally=dev.off())
  }

  mix_p_values <- transform(subset(p_values, is_control & !is_wt_control & include_mix == 1 & perc_expect < 100.0),
                            label = sprintf('%s mutant, genotype %d', mixture_perc_expected, geno_num))
  # This is so to control axis order
  mix_p_values <- transform(mix_p_values, label = factor(label, levels=rev(levels(label))))

  write.csv(mix_p_values, '2_prep_table_mixtures.csv', row.names=FALSE)
  d <- subset(mix_p_values, location %in% c('K65R', 'M184VI'))
  p(d, '2a_prep_table_mixtures_exppos')
  pdf('2c_prep_table_k65r_m184vi.pdf', width=output_width * 2/3, height=output_width * 2/3)
  tryCatch(print(mixture_level_plot(d)), finally=dev.off())
  d <- subset(mix_p_values, location %in% c('K65R', 'M184VI', 'K103NSTH'))
  p(d, '2b_prep_table_mixtures_exppos_k103')


  pdf('prep_table_mixtures_other.pdf', width=output_width * 1.3, height=output_width)
  tryCatch({
    d <- subset(mix_p_values, location %in% c('V106AM', 'Y188LCH', 'G190SAEQ', 'T215YF',
                                              'D67N', 'K101EPH', 'K70E', 'V179DEF'))
    print(mixture_plot(d))
    print(mixture_singlepos_plot(d))
  }, finally=dev.off())
})


#################
# Wild-type plots
#################
local({
    d <- subset(p_values, !is.na(wt_single_barcode) & wt_single_barcode & location %in% prep_primary_sites)
    d <- transform(d, group=sub(' [^ ]+$', '', sample_name))
    d <- transform(d, group=sub('^2059wt', 'MB2059_pol', group))
    d <- transform(d, group=sub('(_wt|-\\d)$', '', group))
    d <- arrange(d, group, subj_id)
    d_ply(d, .(group), function(piece) {
      pdf(sprintf('prep_table_wt_%s.pdf', piece$group[1]), width=output_width * 0.6, height=output_width * 0.03 * length(unique(piece$subj_id)))
      tryCatch(print(wt_plot(piece)), finally=dev.off())
    })

    pdf('7_prep_table_wt_error_by_site_genotype.pdf', width=8, height=5)
    tryCatch({
      print(wt_error_by_site(d))
    }, finally=dev.off())

    # Sum across
    d <- subset(p_values, !is.na(wt_grouped) & wt_grouped & location %in% prep_primary_sites)
    sample <- gsub('(^\\d-|_?wt|_pol|-\\d$|MB)', '', d$sample_name)
    d$genotype <- factor(sample,
                         levels = paste('2059', c('', '_G', '_H', '_E', '_F'), sep=''),
                         labels = c('1', '2', '3', '4', '5'))

    write.csv(wt_table(d),
              '1a_prep_table_wt_avg.csv',
              row.names=FALSE)
    write.csv(d, '1a_prep_table_wt_src.csv', row.names=FALSE)
    write.csv(wt_table(d, function(x) sprintf("%.2f%% - %.2f%%", min(x * 100), max(x * 100))),
              '1b_prep_table_wt_range.csv',
              row.names=FALSE)

    pdf('8_prep_table_wt_grouped_false_positive_rate.pdf', width=output_width, height=length(unique(d$subj_id)) * 15/50)
    tryCatch({
      print(sample_plot(d) + ylab('Control') +
            facet_grid(genotype ~ location, scales = "free_y", space = "free") +
            xlab(sprintf("False positive rate: %.2f%%",
                         100 * sum(d$p_adj < 0.05) / nrow(d))))
      print(sample_pvalue_plot(d) + ylab('Control') +
            facet_grid(genotype ~ location, scales = "free_y", space = "free") +
            xlab(sprintf("False positive rate: %.2f%%",
                         100 * sum(d$p_adj < 0.05) / nrow(d))))
    }, finally=dev.off())
})


#############################################################
# By subject - restricted to subjects with 1+ positive sample
#############################################################

# For each subject, determine whether any loci are positive
any_positive <- ddply(to_plot, .(subj_id), summarize,
                      any_positive=any(!is.na(above_detection) & above_detection))
lu <- match(to_plot[, 'subj_id'], any_positive[, 'subj_id'])
pos <- to_plot[any_positive[lu, 'any_positive'], ]

# Re-calculate
local({
  n_subjects <- length(unique(pos$subj_id))
  write.csv(pos, '5a_prep_table_pos.csv', row.names=FALSE)
  ggsave('5a_prep_table_pos.pdf', sample_plot(pos), width=output_width, height=n_subjects * 15/50)
  ggsave('5b_prep_table_pos_p.pdf', sample_pvalue_plot(pos), width=output_width, height=n_subjects * 15/50)
})

# Splitting M184VI
local({
  l <- c('M184VI', 'M184I', 'M184V')
  to_plot <- subset(p_values, !is_control & visit %in% c('SC', 'SC+1') & location %in% l)
  to_plot$location <- factor(to_plot$location, levels=l)
  to_plot <- annotate_missing_visits(to_plot)
  lu <- match(to_plot[, 'subj_id'], any_positive[, 'subj_id'])
  pos <- to_plot[any_positive[lu, 'any_positive'], ]
  n_subjects <- length(unique(pos$subj_id))
  output_width <- output_width + 2
  write.csv(pos, '6_prep_table_pos_vi.csv', row.names=FALSE)
  ggsave('6_prep_table_pos_vi.pdf', sample_plot(pos), width=output_width, height=n_subjects * 15/50)
  ggsave('6_prep_table_pos_vi_p.pdf', sample_pvalue_plot(pos), width=output_width, height=n_subjects * 15/50)
})

############################################
# FDR
############################################
local({
  l <- c('M184V', 'K65R', 'Y181CIVSFG', 'K103NSTH')
  to_plot <- subset(p_values,
                    ((!is_wt_control & is_control) |
                     (!is.na(wt_grouped) & wt_grouped)) &
                    location %in% l)
  to_plot <- transform(to_plot, type = ordered(ifelse(is_wt_control, 'wt', 'mixture'), levels=c('wt', 'mixture')),
                       called_positive = p_adj < 0.05)
  (xt <- xtabs(~type + called_positive + location, to_plot, drop.unused.levels=TRUE))

  # FDR
  r <- ldply(1:length(l), function(i) {
    loc <- dimnames(xt)$location[i]
    v <- as.vector(xt[, , i])
    data.frame(location=loc, TN=v[1], FN=v[2], FP=v[3], TP=v[4])
  })
  r$FDR <- r$FP / (r$FP + r$TP)
  cs <- colSums(r[, -1])
  overall_fdr <- cs['FP'] / (cs['FP'] + cs['TP'])
})


############################################
# Finally, write a map of the obfuscated IDs
############################################
write.csv(to_plot, 'prep_table_lookup.csv', row.names=FALSE)
