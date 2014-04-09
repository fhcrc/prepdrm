#!/usr/bin/env Rscript

infiles <- commandArgs(TRUE)

raw_frames <- lapply(as.list(infiles), read.csv, as.is=TRUE)

# Find union of all column names
all_headers <- Reduce(function(x, y) c(x, setdiff(y, x)),
                      Map(colnames, raw_frames))

# Add empty columns to files missing some subset
frames <- lapply(raw_frames, function(df) {
  headers <- colnames(df)
  df[,setdiff(all_headers, headers)] <- NA
  df
})

result <- do.call('rbind', frames)

write.csv(result, stdout(), row.names=FALSE, na='')
