#!/usr/bin/env Rscript

library(plyr)
library(ggplot2)
library(RColorBrewer)

# Primary outcomes
prep_primary_sites <- c('K65R', 'K70E', 'M184VI')
# Sites with mutations in controls
prep_control_sites <- c('K65R', 'K103NS', 'M184VI', 'Y181CIV')
# Pattern to match sites which use a genotype-specific correction
genotype_specific_sites <- '^(K65|D67)'

# Global Style
output_width <- 8.5
theme_set(theme_bw(16) + theme(panel.grid = element_blank(),
                               strip.background = element_blank()))

# Below detection fill color
bd_color <- '#AAAAAA'
# Missing data fill color
na_color <- '#AAAAAA'

# Palette for filling cells
fill_colors <- c(bd_color, brewer.pal(8, 'Blues')[-(1:3)])
fill_scale <- scale_fill_manual("% mutant", values=fill_colors, na.value=na_color, drop=FALSE)

# Text colors - dark and light, depending on fill
text_colors <- c('#333333', '#EEEEEE')
text_scale <- scale_color_manual(values=text_colors, guide=FALSE, na.value='#EEEEEE')

# Palette for filling cells colored by p-value
p_value_fill <- scale_fill_manual('p-value\n(FDR corrected)',
                                  values=rev(c(bd_color, brewer.pal(4, 'OrRd')[-1])),
                                  drop=FALSE,
                                  na.value=na_color)


categorize_p_values <- function(p) {
  cut(p,
      c(0, 0.001, 0.01, 0.05, 1),
      labels=c('< 0.001', '< 0.01', '< 0.05', '> 0.1'),
      include.lowest=TRUE)
}

# Plotting functions
sample_plot <- function(df) {
  ggplot(df, aes(x=visit, y=subj_id, fill=mut_prop_cat)) +
    geom_tile() +
    geom_text(aes(label=mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, mut_prop >= 0.1))) +
    facet_wrap(~location) +
    fill_scale +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Subject') +
    xlab('Visit')
}

sample_pvalue_plot <- function(df) {
  df$p_value_cat <- categorize_p_values(df[, 'p_adj'])
  df$p_value_cat[is.na(df$mut_prop_cat)] <- NA
  ggplot(df, aes(x=visit, y=subj_id, fill=p_value_cat)) +
    geom_tile() +
    geom_text(aes(label=mut_prop_label, color=ifelse(is.na(mut_prop_cat), NA, p_adj < 0.01))) +
    facet_wrap(~location) +
    p_value_fill +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Subject') +
    xlab('Visit')
}

mixture_plot <- function(df) {
  df$p_value_cat <- categorize_p_values(df[, 'p_adj'])
  ggplot(df, aes(x = location, y = sample_name, fill = p_value_cat)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100), color = p_adj < 0.01)) +
    p_value_fill +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Mixture') +
    xlab('Location')
}

mixture_singlepos_plot <- function(df) {
  ggplot(df, aes(x = location, y = sample_name, fill = p_adj < 0.05 & n_with_mut == 1)) +
    geom_tile() +
    geom_text(aes(label = sprintf('%.2f%%', mut_prop * 100), color = p_adj < 0.05 & n_with_mut == 1)) +
    scale_fill_manual(values=c('white', brewer.pal(3, 'Blues')[3]), guide=FALSE) +
    text_scale +
    theme(panel.grid = element_blank(),
          strip.background = element_blank()) +
    ylab('Mixture') +
    xlab('Location')
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


########################
## Load and plot data ##
########################
p_values <- read.csv('output/prep_merged_fisher.csv', as.is=TRUE)

## Generate an anonymized subject id
subj_id_lu <- ddply(p_values, .(sample_name), summarize,
                    subj_id = factor(sub('p(\\d)d.bc(\\d{3}).*', '\\1-\\2', sample_ids[1])))
subj_id_lu <- mutate(subj_id_lu, subj_id = factor(subj_id, levels=rev(sort(levels(subj_id)))))

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
                   # Assign an ID: plate and first barcode sample
                   subj_id = subj_id_lu$subj_id[match(sample_name, subj_id_lu$sample_name)],
                   # Recode SC+3 as SC+1, recode SC+1repeat as SC+1
                   visit = sub(' ?repeat', '', ifelse(visit == 'SC+3', 'SC+1', visit)),
                   # Order locations by position
                   location = reorder(location, as.integer(sub('^[A-Z]+(\\d+).*$', '\\1', as.character(location))), max))

# Add bd category
p_values <- mutate(p_values,
                   mut_prop_cat = ordered(ifelse(!above_detection, 'bd',
                                                 as.character(mut_prop_cat)),
                                          levels=c('bd', levels(mut_prop_cat))))

# Subset to samples, sites of interest
to_plot <- subset(p_values, !is_control & visit %in% c('SC', 'SC+1') & location %in% prep_primary_sites)
# Add rows with missing values for missing samples
to_plot <- ddply(to_plot, .(location, subj_id), function(piece) {
  for(v in c('SC', 'SC+1')) {
    if(!(v %in% piece$visit)) {
      # Add row
      piece <- rbind(piece, mutate(piece[1, ], visit=v, mut_prop_label='no data', mut_prop_cat=NA, above_detection=NA))
    }
  }
  return(piece)
})

# Equivalent for mixture controls
wt_to_plot <- subset(p_values, is_wt_control & location %in% prep_primary_sites)

# Number of subjects
n_subjects <- length(unique(to_plot$subj_id))

######################
# Big by subject table
######################
pdf('prep_table_all.pdf', width=output_width, n_subjects * 12/50)
tryCatch(print(sample_plot(to_plot)), finally=dev.off())

n_splits <- 2
pdf('prep_table_all_split.pdf', width=output_width, n_subjects * 12/50 / n_splits)
tryCatch({
  pv <- mutate(to_plot, split=cut(as.integer(to_plot$subj_id), breaks=n_splits))
  d_ply(pv, .(split), function(piece) {
    print(sample_plot(piece))
  })
}, finally=dev.off())


#################
# Mixture plots
#################
pdf('prep_table_mixtures_exppos.pdf', width=output_width, height=output_width)
tryCatch({
  d <- subset(p_values, is_control & !is_wt_control & location %in% prep_control_sites)
  print(mixture_plot(d))
  print(mixture_singlepos_plot(d))
}, finally=dev.off())
pdf('prep_table_mixtures_other.pdf', width=output_width * 1.3, height=output_width)
tryCatch({
  d <- subset(p_values, is_control & !is_wt_control & location %in% c('V106AM', 'Y188LCH', 'G190SAEQ', 'T215YF',
                                                                      'D67N', 'K101EPH', 'K70E', 'V179DEF'))
  print(mixture_plot(d))
  print(mixture_singlepos_plot(d))
}, finally=dev.off())


#################
# Wild-type plots
#################
local({
    d <- subset(p_values, is_wt_control & location %in% prep_primary_sites)
    d <- transform(d, group=sub(' [^ ]+$', '', sample_name))
    d <- transform(d, group=sub('^2059wt', 'MB2059_pol', group))
    d <- transform(d, group=sub('_wt$', '', group))
    d <- arrange(d, group, subj_id)
    d_ply(d, .(group), function(piece) {
      pdf(sprintf('prep_table_wt_%s.pdf', piece$group[1]), width=output_width * 0.6, height=output_width * 0.03 * length(unique(piece$subj_id)))
      tryCatch(print(wt_plot(piece)), finally=dev.off())
    })
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
n_subjects <- length(unique(pos$subj_id))
ggsave('prep_table_pos.pdf', sample_plot(pos), width=output_width, height=n_subjects * 15/50)
ggsave('prep_table_pos_p.pdf', sample_pvalue_plot(pos), width=output_width, height=n_subjects * 15/50)

############################################
# Finally, write a map of the obfuscated IDs
############################################
write.csv(to_plot, 'prep_table_lookup.csv', row.names=FALSE)
