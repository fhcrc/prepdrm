#!/usr/bin/env Rscript

closest <- read.csv('../docs/closest_reference.csv', as.is=TRUE)

files <- sort(list.files(pattern='plate\\d+.csv'))

for(f in files) {
  message(f)
  df <- read.csv(f, as.is=TRUE, comment.char='#')
  mc <- closest$min_dist[match(df$reference_id, closest$sequence)]
  df$matched_control <- ifelse(df$matched_control == '', mc, df$matched_control)
  write.csv(df, f, row.names=FALSE, na='')
}
