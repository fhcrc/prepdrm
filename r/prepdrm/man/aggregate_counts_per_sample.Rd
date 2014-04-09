\name{aggregate_counts_per_sample}
\alias{aggregate_counts_per_sample}
\title{Aggregate counts per sample_name visit combination}
\usage{
aggregate_counts_per_sample(counts, grouping = .(location, sample_name, visit,
  matched_control, is_control, is_wt_control))
}
\arguments{
  \item{counts}{Data frame, as from
  \code{load_counts_and_meta}}

  \item{grouping}{Grouping values}
}
\value{
Data frame
}
\description{
Aggregate counts per sample_name visit combination
}

