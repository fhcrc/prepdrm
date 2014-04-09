#!/usr/bin/env Rscript
library(methods)
library(devtools)
library(plyr)
library(reshape2)
load_all("r/prepdrm", export_all = FALSE)

args <- commandArgs(TRUE)

if (length(args) != 4) {
  stop("USAGE: aggregate_aa_per_sample.R meta count output output_wide")
}

r <- aggregate_aa_per_sample(load_aa_and_meta(args[2], args[1]))


write.csv(r, args[3], row.names = FALSE)

r_long <- ddply(r, .(location, sample_name, visit), function(piece) {
  loc <- strsplit(piece$location[1], c())[[1]]
  m <- melt(piece, id.vars = 1:3, variable.name = "amino_acid")
  m[as.character(m$amino_acid) %in% loc, ]
})

r_wide <- dcast(fix_visits(r_long), sample_name ~ location + amino_acid + visit, fill=0)

write.csv(r_wide, args[4], row.names = FALSE)
