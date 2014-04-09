#!/usr/bin/env Rscript
library(reshape2)

args <- commandArgs(TRUE)

if(length(args) != 3)
  stop("USAGE: sequence_classification_by_sample.R infile metadata outfile")
input <- read.csv(args[1], as.is=TRUE)
meta <- read.csv(args[2], as.is=TRUE)
colnames(meta)[1] <- 'sample'

m <- merge(input, meta[, c('sample', 'subject', 'visit', 'direction')], all.x=TRUE)

melted <- melt(m, measure.vars='count')
casted <- dcast(melted, subject+visit+location+amino_acid+classification+
                        sequence+codon_highlighted+any_indels+whole_codon_indels~variable+direction,
                sum, fill=0)

casted <- with(casted, casted[order(subject, visit, location, -count_forward-count_reverse),])

write.csv(casted, args[3], row.names=FALSE, quote=FALSE, na='')
