\name{sw.plot}
\alias{sw.plot}
\title{Plot Results of Smith-Waterman Algorithm}
\description{
  This function plots the sign-adjusted logratios by their chromosomal location.
    It can superimpose the location of the highest-scoring island found by the
    Smith-Waterman algorithm, the results of a robustness analysis, and the
    expected logratios based on known copy numbers in the test DNA.
  }
\usage{
  sw.plot(logratio, location = seq(length(logratio)),
    threshold.func = function(x) median(x) + .2 * mad(x),
    sign = -1, highest = TRUE, expected = NULL, rob = NULL, legend = TRUE, 
    xlab = "Chromosomal location", ylab = "Intensity log ratio", ...)
  }
\arguments{
  \item{logratio}{a vector of logratios, not adjusted for sign or threshold}
  \item{location}{a vector of chromosomal locations corresponding to the log ratios}
  \item{threshold.func}{threshold function: see \code{\link{sw.threshold}}}
  \item{sign}{sign of logratio adjustment: see \code{\link{sw.threshold}}}
  \item{highest}{plot location of highest-scoring island if TRUE}
  \item{expected}{a vector of expected copy numbers, or NULL}
  \item{rob}{a vector of robustness scores, or NULL}
  \item{legend}{plot legend if TRUE}
  \item{xlab}{X axis label}
  \item{ylab}{Y axis label}
  \item{\dots}{other arguments passed to the 'plot' function}
  }
\author{T.S.Price}
\seealso{
  \code{\link{sw}}
  \code{\link{sw.threshold}}
  \code{\link{sw.perm.test}}
  \code{\link{sw.rob}}
  }
\examples{
## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## perform permuation test for islands identified
p <- sw.perm.test(x, max.nIslands = NULL, nIter = 1e4)

## calculate robustness scores
r <- sw.rob(x)

## plot results
sw.plot(logratio, seq(length(logratio)),
  function(x) median(x) + .2 * mad(x), sign = -1, rob = r,
  main = paste("Toy dataset, highest-scoring island p =", p[1]))
  }
\keyword{misc}
