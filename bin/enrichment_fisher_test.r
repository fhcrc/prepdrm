#!/usr/bin/env Rscript

# Test for mutation enrichment using Fisher's exact test

library(ggplot2)
library(plyr)

library(methods)
library(devtools)
load_all('r/prepdrm', export_all=FALSE)

theme_set(theme_bw(16))

args <- commandArgs(TRUE)
stopifnot(length(args) == 5)

sample_meta_path <- args[1]
counts_path <- args[2]
long_output_path <- args[3]
wide_output_path <- args[4]
pdf_output_path <- args[5]

message('loading counts and metadata')
m <- load_counts_and_meta(counts_path, sample_meta_path)

## Only use matched controls for K65, D67
genotype_spec <- grepl(prepdrm::genotype_specific_sites, m$location)
m$reference_id[m$is_wt_control & !genotype_spec] <- ''
m$matched_control[!genotype_spec] <- ''

## Control counts
wt_controls <- subset(m, is_wt_control)
controls <- ddply(wt_controls, .(reference_id, location),
                  summarize,
                  wt=sum(wt), mut=sum(mut), other=sum(other),
                  non_coding=sum(non_coding), nonmutant=sum(nonmutant), sample_ids=paste(sample_id, collapse=','))
colnames(controls)[-(1:2)] <- paste('control', colnames(controls)[-(1:2)], sep='_')

message('aggregating counts per-sample')
cases <- aggregate_counts_per_sample(subset(m, !is_wt_control))

message('Match controls with cases')
cases_controls <- merge(cases, controls,
                        by.y=c('location', 'reference_id'),
                        by.x=c('location', 'matched_control'), all.x=TRUE)

each_wt_control <- ddply(wt_controls, .(reference_id, location), function(piece) {
    ddply(piece, .(sample_id), function(barcode) {
      matched <- summarize(piece[piece$sample_id != barcode$sample_id, ],
                           wt=sum(wt), mut=sum(mut), other=sum(other),
                           non_coding=sum(non_coding), nonmutant=sum(nonmutant), sample_ids=paste(sample_id, collapse=','))
      colnames(matched) <- paste('control', colnames(matched), sep='_')
      stopifnot(nrow(matched) == 1)
      stopifnot(nrow(barcode) == 1)
      m <- cbind(barcode, matched)
      stopifnot(nrow(m) == 1)
      transform(m, n_with_mut=1, prop_with_mut=1, mut_counts='', sample_ids=sample_id,
                sample_name=paste(sample_name, sample_id))
    })
})

each_wt_control <- each_wt_control[, setdiff(colnames(each_wt_control),
                                             c('reference_id', 'sample_id',
                                               'plate', 'barcode','direction',
                                               'sample', 'wt_prop',
                                               'other_prop',
                                               'non_coding_prop'))]
each_wt_control$wt_single_barcode <- TRUE

# Test groups of 4 wt controls against other controls
grouped_wt_control <- ddply(aggregate_counts_per_sample(transform(wt_controls,
                                                                  sample_name=sprintf('%d-%s', plate, sample_name))),
                            .(reference_id, location),
                            function(piece) {
    ddply(piece, .(sample_name), function(sample) {
      matched <- summarize(piece[piece$sample_name != sample$sample_name, ],
                           wt=sum(wt), mut=sum(mut), other=sum(other),
                           non_coding=sum(non_coding), nonmutant=sum(nonmutant), sample_ids=paste(sample_ids, collapse=','))
      colnames(matched) <- paste('control', colnames(matched), sep='_')
      stopifnot(nrow(matched) == 1)
      stopifnot(nrow(sample) == 1)
      m <- cbind(sample, matched)
      stopifnot(nrow(m) == 1)
      m
    })
})

grouped_wt_control <- grouped_wt_control[, setdiff(colnames(grouped_wt_control),
                                                   c('reference_id', 'sample_id',
                                                     'plate', 'barcode','direction',
                                                     'sample', 'wt_prop',
                                                     'other_prop',
                                                     'non_coding_prop'))]
grouped_wt_control$wt_grouped <- TRUE

cases_controls <- rbind.fill(cases_controls, each_wt_control, grouped_wt_control)

samples_controls <- merge(m[!m$is_wt_control, ], controls,
                          by.y=c('location', 'reference_id'),
                          by.x=c('location', 'matched_control'), all.x=TRUE)

# Run test
test_elevated <- function(control_mut, control_nonmut, mut, nonmut) {
  tab <- matrix(c(control_nonmut, control_mut, nonmut, mut), byrow=TRUE, nrow=2,
                dimnames=list(c('control', 'sample'), c('nonmutant', 'mutant')))
  fisher.test(tab, alternative='greater', conf.int=TRUE)
}

message('Testing for mutation enrichment')
sample_elevated <- function(piece) {
  if(nrow(piece) != 1)
    print(piece)
  stopifnot(nrow(piece) == 1)
  piece <- transform(piece, mut_prop=mut/(mut+nonmutant),
                     control_mut_prop=control_mut / (control_mut+control_nonmutant))
  if(is.na(piece$control_mut))
    return(transform(piece,
                     p_value=NA,
                     or=NA,
                     alternative=NA))

  r <- with(piece, test_elevated(control_mut, control_nonmutant, mut, nonmutant))
  transform(piece,
            p_value=r$p.value,
            or=r$estimate,
            alternative=r$alternative)
}
sample_p_values <- ddply(samples_controls,
                         .(location, sample_id), sample_elevated,
                         .progress='text')
message('Sample positive rates')
sample_pos <- ddply(sample_p_values, .(location, sample_name, visit),
                    function(piece) {
  data.frame(n_pos_0.05 = sum(piece$p_value < 0.05),
             p_values = paste(piece$p_value, collapse=', '))
}, .progress='text')

p_values <- ddply(cases_controls, .(location, sample_name, visit), sample_elevated)
p_values <- merge(p_values, sample_pos, by=c('location', 'sample_name', 'visit'),
                  all.x=TRUE)

# Adjusted for control samples only, only for primary mutations
p_values$p_adj <- adjust_p_values(p_values, controls_only = FALSE)
write.csv(p_values, long_output_path, row.names=FALSE, na='')

# Wide
message('Wide results')
write_wide_results(p_values, wide_output_path, extra_vars='p_adj')

message('Plotting')
warning('Skipping plots for control only set')
#write_enrichment_plot(p_values, pdf_output_path)
