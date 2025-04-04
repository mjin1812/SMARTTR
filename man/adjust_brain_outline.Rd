% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping_functions.R
\name{adjust_brain_outline}
\alias{adjust_brain_outline}
\title{Adjust brain outline.}
\usage{
adjust_brain_outline(s, filter = NULL)
}
\arguments{
\item{s}{slice object}

\item{filter}{(list, default = NULL) If the user passes their own filter list, it will use that instead of the presaved filter list from SMARTTR.}
}
\value{
filter (list) wholebrain compatible filter
}
\description{
This function takes a slice object and first applies a filter with default settings to
the image set as the slice registration path. An interative user loop allows for easy adjustment of the brain threshold
since the wholebrain GUI tends to be a bit buggy. This function then returns a filter with the adjusted brain threshold.
}
\examples{
\dontrun{
# Adjust the brain threshold, then run register on the slice object
filter <- adjust_brain_outline(s); s <- register(s, filter = filter)
}
}
