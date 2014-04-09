#!/usr/bin/env Rscript

library(ggplot2)
library(plyr)
library(reshape2)

theme_set(theme_bw(16))

args <- commandArgs(TRUE)
stopifnot(length(args) == 2)

infile <- args[1]
outfile <- args[2]

pat <- 'p(\\d+)d([12])bc(\\d+)'
err_rate <- transform(read.csv(infile, as.is=TRUE),
                      plate=sub(pat, '\\1', sample),
                      direction=sub(pat, '\\2', sample))
err_rate <- subset(err_rate, direction=='1')

# Drop samples without reference
err_rate <- subset(err_rate, !grepl('p2d.bc340|p3d.bc317', sample))

m <- melt(err_rate, measure.vars=c('n_del', 'n_ins', 'n_match', 'n_mismatch'))
m <- transform(m, proportion=value/n_aligned,
               type=ifelse(grepl('n_(del|ins)', as.character(variable)), 'Insertion / Deletion', 'Aligned Base'))
m <- transform(m, variable=sub('n_', '', as.character(variable)))

pdf(outfile, width=10, height=6)
d_ply(m, .(feature, ref_sequence), function(piece) {
  labels <- unique(piece[, c('position', 'ref_base')])
  labels <- labels[order(labels$position),]
  n_reads <- with(subset(piece, position==min(position) & variable == 'match'), sum(n_aligned))
  message(n_reads)
  p <- ggplot(piece, aes(x=ordered(position+1), y=proportion, fill=variable)) +
    geom_boxplot() +
    ggtitle(paste(piece$feature[1], piece$ref_sequence[1], sep='\n')) +
    scale_x_discrete(name="Reference Base", breaks=ordered(labels$position+1), labels$ref_base) +
    ylab("Proportion of reads") +
    scale_fill_discrete(name="Type") +
    facet_wrap(~type, ncol=1)
  print(p)
})
dev.off()
