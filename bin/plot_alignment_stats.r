#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)

theme_set(theme_bw(16))

args <- commandArgs(TRUE)
if(length(args) != 3)
  stop('usage: plot_alignment_stats.R stats meta out_pdf')

infile <- args[1]
meta_path <- args[2]
outfile <- args[3]

d <- read.csv(infile, as.is=TRUE)
meta <- read.csv(meta_path, as.is=TRUE)
d <- merge(d, meta, by.x=1, by.y='sample_id', all.x=TRUE)

variable_names <- data.frame(name=c('mean_aligned_length', 'mean_pct_id', 'n'),
                             label=c('Mean aligned sequence length', 'Mean % ID', '# of reads'),
                             stringsAsFactors=FALSE)

m <- melt(d, id.vars=c('sample', 'plate', 'barcode_id', 'is_control', 'direction'),
          measure.vars=variable_names$name)
m <- transform(m, label=variable_names$label[match(variable, variable_names$name)])

p1 <- ggplot(m, aes(x=factor(plate), y=value, fill=factor(direction))) +
  geom_boxplot() +
  facet_wrap(~label, scales='free_y') +
  theme(legend.position='bottom') +
  xlab('Run') +
  ggtitle('Alignment statistics by run')


p2_d <- subset(d, plate > 1) #subset(d, !is_control)
n_reads_by_barcode <- aggregate(n~barcode_id, p2_d, median)
n_reads_by_barcode <- n_reads_by_barcode[order(n_reads_by_barcode$n),]
p2 <- ggplot(p2_d, aes(x=factor(barcode_id, levels=n_reads_by_barcode$barcode_id), y=n)) +
  geom_boxplot() +
  geom_point(aes(color=factor(plate))) +
  theme_bw(10) +
  theme(axis.text.x=element_text(angle=90)) +
  theme(legend.position='bottom') +
  xlab('Barcode #') +
  ylab('# of reads')

pdf(outfile, width=14, height=5)
print(p1)
print(p2)
dev.off()
