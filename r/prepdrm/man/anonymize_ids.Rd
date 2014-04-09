\name{anonymize_ids}
\alias{anonymize_ids}
\title{Assign a unique subject ID using first plate, barcode for each subject}
\usage{
anonymize_ids(sample_ids, sample_names)
}
\arguments{
  \item{sample_names}{Sample names - for grouping}

  \item{sample_ids}{Sample IDs, i.e. p1d1bc241}
}
\value{
a new sample_ids
}
\description{
Assign a unique subject ID using first plate, barcode for
each subject
}

