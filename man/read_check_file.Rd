% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{read_check_file}
\alias{read_check_file}
\title{Read a csv or excel file as a tibble. Checks first that the file exists, and that it is a csv or xlsx format.}
\usage{
read_check_file(x, ...)
}
\arguments{
\item{x}{A file path}

\item{...}{additional parameters to pass to either the \link[readr:read_delim]{readr::read_csv} function or \link[readxl:read_excel]{readxl::read_excel}}
}
\value{
tibble dataframe
}
\description{
Read a csv or excel file as a tibble. Checks first that the file exists, and that it is a csv or xlsx format.
}
