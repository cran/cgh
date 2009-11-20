\name{sw.perm.test}
\alias{sw.perm.test}
\title{Permutation Test for Smith-Waterman Algorithm}
\description{Perform a permutation test of island scores from the Smith-Waterman algorithm.}
\usage{
  sw.perm.test(x, max.nIslands = 1, nIter = 1000, seed = NULL, trace = FALSE)
  }
\arguments{
  \item{x}{a vector of real values}
  \item{max.nIslands}{number of iterations of the algorithm,
    each iteration finding the next highest-scoring 'island' of
    positive values, or NULL to find all islands}
  \item{nIter}{number of permutations of the input data used in the test}
  \item{seed}{seed for the random number generator, or NULL to use a faster
    random number generator that cannot be seeded}
  \item{trace}{print verbose output if TRUE}
  }
\value{A vector of probability values, calculated as the proportion of instances for which
  performing the Smith-Waterman algorithm on random permutations of the data
  identifies a higher-scoring island than the islands identified when the algorithm
  is performed on the original data}
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
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## perform Smith-Waterman
sw(x)

## perform permutation test on the islands identified
sw.perm.test(x, max.nIslands = NULL, nIter= 1e4)
  }
\keyword{misc}
