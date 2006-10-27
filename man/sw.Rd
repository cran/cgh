\name{sw}
\alias{sw}
\title{Perform the Smith-Waterman Algorithm}
\description{Perform the Smith-Waterman algorithm on a vector of real values.}
\usage{
  sw(x, max.nIslands = NULL, trace = FALSE)
  }
\arguments{
  \item{x}{a vector of real values}
  \item{max.nIslands}{the number of iterations of the algorithm performed.
    Each iteration finds the next highest-scoring 'island' of
    positive values. Set to NULL to find all islands}
  \item{trace}{print verbose output if TRUE}
  }
\details{
  The Smith-Waterman algorithm detects 'islands' of positive scores in a vector of
    real values. The input values should have a negative mean. The algorithm can be
    used to identify regions of copy number change in microarray fluorescence
    logratios, once the logratios have been adjusted for sign and a suitable
    threshold value subtracted to ensure a negative mean: see \code{\link{sw.threshold}}}
\value{
  \item{x}{the input vector}
  \item{s}{a numeric vector containing the partial sums after one iteration
    of the Smith-Waterman algorithm}
  \item{score}{a numeric vector of island scores}
  \item{start}{a numeric vector of indices identifying the start of each island}
  \item{length}{a numeric vector of island lengths}
  }
\references{
  Smith TF, Waterman MS.
  Identification of common molecular subsequences.
  J Mol Biol. 1981;147(1):195-7.
  }
\author{T.S.Price}
\seealso{
  \code{\link{sw.threshold}}
  \code{\link{sw.perm.test}}
  \code{\link{sw.rob}}
  \code{\link{sw.plot}}
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
