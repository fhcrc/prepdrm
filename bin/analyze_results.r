#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
library(plyr)

theme_set(theme_bw(14))

args <- commandArgs(TRUE)

stopifnot(length(args) == 3)
infile <- args[1]
outfile <- args[2]
out_csv <- args[3]

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

mutation_location_factor <- function(locations) {
  l_names <- unique(locations)
  locs <- as.integer(sub('^[A-Z]+(\\d+).*$', '\\1', l_names))
  l_names <- l_names[order(locs)]
  return(ordered(locations, levels=l_names))
}

pat <- 'p(\\d+)d([12])bc(\\d+)'

mut <- read.csv(infile, as.is=TRUE)
mut <- transform(mut,
                 barcode=ordered(as.integer(sub(pat, '\\3', sample))),
                 plate=ordered(as.integer(sub(pat, '\\1', sample))),
                 #plate_barcode=factor(sub(pat, '\\1-\\3', sample)),
                 direction=ifelse(grepl('p\\d+d2', sample), 'reverse', 'forward'),
                 classification=ordered(classification, levels=c('wt', 'mut', 'other', 'non_coding')),
                 location=mutation_location_factor(location))
agg <- aggregate(count~location+sample+plate+barcode+direction+classification, mut, sum)
read_counts <- aggregate(count~location+sample+plate+barcode+direction, mut, sum)

scale <- scale_fill_manual(name='Class', values=c(gg_color_hue(3), 'grey'))

dfs <- list(full=agg)

#dfs <- list(full=agg, subset=subset(agg, classification != 'non_coding'))
#stopifnot(length(dfs) == 2)

# 12/84 inches per mutation/barcode
#height <- 12 * length(unique(mut$location))
width <- 12 * nrow(unique(mut[,c('barcode', 'plate')])) / 84

pdf(outfile, width=width, height=12)

# Read counts by location
d_ply(read_counts, .(location), function(loc) {
  p <- ggplot(loc, aes(x=barcode, y=count)) +
    geom_bar(stat='identity') +
    facet_grid(direction~plate, scales='free_x', space='free_x') +
    theme(axis.text.x=element_text(angle=-90)) +
    ggtitle(paste("Location coverage:", loc$location[1])) +
    xlab("Barcode Number") +
    ylab("Read count")

  print(p)
})

# Mutation proportions by location
l_ply(dfs, function(df) {
  d_ply(df, .(location), function(loc) {
    p <- ggplot(loc, aes(x=barcode, y=count, fill=classification)) +
      geom_bar(stat='identity', position='fill') +
      scale +
      theme(axis.text.x=element_text(angle=-90)) +
      ggtitle(paste("Mutation frequencies:", loc$location[1])) +
      xlab("Barcode Number") +
      ylab("Proportion") +
      theme(legend.position='bottom') +
      facet_grid(direction~plate, scales='free_x', space='free_x')
    print(p)
  }, .progress='text')
})
dev.off()

# Table
m <- melt(agg, measure.vars='count')
cst <- dcast(m, location+plate+barcode+direction+sample~classification, sum)

# Add props
pt <- prop.table(as.matrix(cst[, -(1:5)]), margin=1)
colnames(pt) <- paste(colnames(pt), 'prop', sep='_')
cst <- cbind(cst, pt)
write.csv(cst, out_csv, row.names=FALSE)
