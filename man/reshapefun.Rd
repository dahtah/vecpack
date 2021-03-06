% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpack.R
\name{reshapefun}
\alias{reshapefun}
\title{S3 generic for converting an object back into a vector}
\usage{
reshapefun(x, ...)
}
\arguments{
\item{x}{an object}

\item{...}{ignored}
}
\description{
When converting an object into a vector part of the meta-data is lost (for example, when converting a matrix to a vector, the dimensions of the original matrix are lost). These generic functions return a "reshaping function" appropriate for turning a vector back into the original object. They're used internally by vecpack, see Vignette.
}

