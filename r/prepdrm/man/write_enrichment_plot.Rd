\name{write_enrichment_plot}
\alias{write_enrichment_plot}
\title{Visualization of mutation rates at different visits}
\usage{
write_enrichment_plot(df, output_path, locs = c("K65R", "M184VI", "K70E",
  "Y181CIV"), p_value_var = "p_value")
}
\arguments{
  \item{df}{Data frame. Must have columns p_value,}

  \item{output_path}{Output path for PDF}

  \item{locs}{Locations to include}

  \item{p_value_var}{Name of the p value column}
}
\value{
(Invisibly) the plot
}
\description{
Visualization of mutation rates at different visits
}

