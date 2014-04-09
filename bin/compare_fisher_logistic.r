#!/usr/bin/env Rscript

fisher <- read.csv('output/prep_merged_fisher.csv', as.is=TRUE)
logit <- read.csv('output/prep_merged_logistic.csv', as.is=TRUE)

m <- merge(fisher, logit, by=c('location', 'sample_name', 'visit'),
           suffixes=c('_fisher', '_logit'))
