#!/usr/bin/env Rscript

library(methods)
library(devtools)

library(reshape2)
library(plyr)

load_all('r/prepdrm')


usage <- 'usage: prep_subject_results_wide.r (all|primary) meta fisher aa_counts output'
args <- commandArgs(TRUE)


if(length(args) != 5)
  stop(usage)

subset_to_primary <- switch(tolower(args[1]),
                            all=FALSE,
                            primary=TRUE)
if(is.null(subset_to_primary))
  stop(usage)

meta <- read.csv(args[2], as.is=TRUE)
results <- fix_visits(read.csv(args[3], as.is=TRUE))
aa <- read.csv(args[4], as.is=TRUE)

# subset to samples, primary sites
results <- results[!results$is_control, ]
if(subset_to_primary)
 results <- results[results$location %in% prepdrm::prep_primary_sites, ]

# Total reads
results <- mutate(results, total_reads=mut+nonmutant)

# Above 1%
results <- mutate(results,
                  pos_0.05=p_adj < 0.05,
                  pos_0.05_abv1=pos_0.05 & mut_prop >= 0.01)

# Rename a column to make things clearer later
colnames(results) <- sub('n_pos_0.05', 'n_barcode_pos_0.05', colnames(results))

# By visit, across sites
results <- ddply(results, .(sample_name, visit), function(piece) {
  mutate(piece,
         sum_resist_visit=sum(mut_prop[pos_0.05]),
         n_sites_pos_visit=sum(pos_0.05),
         sum_sites_pos_abv1_visit=sum(pos_0.05_abv1),
         any_pos_visit=any(pos_0.05),
         any_pos_abv1_visit=any(pos_0.05_abv1))
})

# Across visits
results <- ddply(results, .(sample_name), function(piece) {
  mutate(piece,
         sum_resist_visit=sum(mut_prop[pos_0.05]),
         n_sites_pos=sum(pos_0.05),
         n_sites_pos_abv1=sum(pos_0.05_abv1),
         any_pos=any(pos_0.05),
         any_pos_abv1=any(pos_0.05_abv1))
})

# Make wide. Do this in two steps: first across sites
r <- results[, grep('sample_name|n_barcode_pos|n_sites|sum_resist|any_pos|^mut_prop|pos_0.05|p_adj|location|visit|total_reads|is_sc3', colnames(results))]
r$matched_control <- meta$matched_control[match(r$sample_name, meta$sample_name)]
wide_by_visit <- reshape(r,
                         direction='wide',
                         timevar='location',
                         idvar=grep('sample_name|visit|matched_control|any_pos|is_sc3|.*_visit|any_pos.*|n_sites.*|sum_resist.*', colnames(r), value=TRUE),
                         sep='_')

# Then across visits
wide_by_subject <- reshape(wide_by_visit,
                           direction='wide',
                           timevar='visit',
                           idvar=grep('^(sample_name|matched_control|any_pos|any_pos_abv1|n_sites_pos|n_sites_posabv1|sum_resist)$', colnames(wide_by_visit),
                                      value=TRUE),
                           sep='_')

# Drop an extra visit label from some of the columns - this was just to prevent a clash earlier
colnames(wide_by_subject) <- sub('_visit_sc', '_sc', colnames(wide_by_subject))
wide_by_subject <- wide_by_subject[, which(colnames(wide_by_subject) != 'is_sc3_sc')]

wide_by_subject <- join(wide_by_subject, aa, by='sample_name', match='first')

# Write the results
write.csv(wide_by_subject, args[5], na='', row.names=FALSE)
