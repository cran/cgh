\name{sw.threshold}
\alias{sw.threshold}
\title{Threshold function}
\description{Function to adjust intensity logratios for sign and threshold before performing the Smith-Waterman Algorithm.}
\usage{
  sw.threshold(logratio,
      threshold.func = function(x) median(x) + .2 * mad(x), sign = +1)
  }
\arguments{
  \item{logratio}{a vector of real values, corresponding to fluorescence intensity logratios}
  \item{threshold.func}{function for calculating threshold}
  \item{sign}{sign of logratio adjustment}
  }
\details{
  The purpose of this function is to adjust the microarray fluorescence intensity logratios to ensure that they have
    the appropriate sign and a mean that is less than zero. sign = +1 is used to detect polysomy (regions of copy number change increase)
    in test:control logratios. Conversely, sign = -1 is used -- inverting the sign of the logratios -- to detect deletions
    (regions of copy number decrease). A threshold, calculated using the threshold function, is subtracted from the sign-adjusted
    logratios to ensure that they have a negative mean. The default threshold function is equal to the median, plus a small
    contant multiplied by a robust estimator of the standard deviation.}
\value{
  A numeric vector equal to

      sign * logratio - threshold.func( sign * logratio )
  }
\author{T.S.Price}
\references{
  Price TS, et al.
  SW-ARRAY: a dynamic programming solution for the identification of copy-number changes in genomic DNA using array comparative genome hybridization data.
  Nucl Acids Res. 2005;33(11):3455-3464.
  }
\seealso{
  \code{\link{sw}}
  }
\examples{
## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), sign = -1)

## perform Smith-Waterman algorithm
sw(x, trace = TRUE)
  }
\keyword{misc}
