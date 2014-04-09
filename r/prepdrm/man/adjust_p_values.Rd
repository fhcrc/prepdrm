\name{adjust_p_values}
\alias{adjust_p_values}
\title{Adjusts p-values}
\usage{
adjust_p_values(df, primary_separate = TRUE, controls_only = TRUE,
  input_column = "p_value", method = "fdr")
}
\arguments{
  \item{df}{Data frame}

  \item{primary_separate}{Correct primary sites
  separately?}

  \item{controls_only}{Only display corrected p-values for
  controls?}

  \item{input_column}{Column containing uncorrected
  p-value}

  \item{method}{passed to \code{p.adjust}}
}
\value{
Corrected p-value
}
\description{
\code{adjust_p_values} adjusts p-values, by default using
FDR correction, for a data frame testing for enriched sites
from PrEP.
}
\details{
By default, only control p-values are shown adjusted,
though the correction uses all. Primary sites (defined in
\code{prep_primary_sites}) are corrected using p-values
only from those sites. Secondary sites are corrected using
all p-values.
}

