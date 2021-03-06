% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/vpack.R
\name{gen.vunpack}
\alias{gen.vunpack}
\title{Generate an unpacking method for vectors}
\usage{
gen.vunpack(l)
}
\arguments{
\item{l}{a list}
}
\value{
an unpacking method
}
\description{
Generate an unpacking method for vectors
}
\examples{
#A list of values that can all be coerced into a vector
l <- list(a= matrix(1:4,2,2),b = 5)
#Pack elements into a vector
vpack(l)
#For the inverse operation (getting the list back from our packed vector)
#we create an unpacking method
unpack <- gen.vunpack(l)
vpack(l) \%>\% unpack
}
\author{
Simon Barthelme
}

