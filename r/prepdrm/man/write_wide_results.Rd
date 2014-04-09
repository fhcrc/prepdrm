\name{write_wide_results}
\alias{write_wide_results}
\title{Convert a long data frame of p-values and results to a wide one, write to a file}
\usage{
write_wide_results(df, output_path, p_value_var = "p_value", cutoff = 0.05,
  extra_vars = c())
}
\arguments{
  \item{df}{Data frame. Must have columns p_value,}

  \item{output_path}{Output path for CSV file}

  \item{p_value_var}{Name of the p value column}

  \item{cutoff}{Cutoff for \code{p_value_var} to consider
  significant}

  \item{extra_vars}{Extra columns to include in teh output}
}
\value{
(Invisibly) the plot
}
\description{
Convert a long data frame of p-values and results to a wide
one, write to a file
}

