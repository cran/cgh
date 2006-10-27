\name{sw.rob}
\alias{sw.rob}
\title{Robustness Calculation for Smith-Waterman Algorithm}
\description{
  Calculate robustness scores to evaluate how sensitive to the threshold value
  is the localisation of the highest-scoring island identified by the Smith-Waterman algorithm
  }
\usage{
  sw.rob(x, lo.func = function(x) median(x),
      hi.func = function(x) median(x) + .4 * mad(x), prec = 100)
  }
\arguments{
  \item{x}{a vector of real values}
  \item{lo.func}{a function for the lowest threshold value}
  \item{hi.func}{a function for the highest threshold value}
  \item{prec}{the precision of the calculation.}
  }
\details{
  This function performs a sensitivity analysis to determine the robustness the localisation
    of the highest-scoring island obtained by the Smith-Waterman algorithm to different values
    of the threshold. The Smith-Waterman algorithm is run repeatedly, each time using a different
    threshold value. The range of threshold values used is that obtained by dividing
    ( lo.func(x), hi.func(x) ) into `prec' equal intervals. The robustness is calculated as the
    proportion of times that a particular chromosomal location falls within the highest-scoring
    island.}
\value{A vector of robustness values equal in length to the input vector.}
\author{T.S.Price}
\seealso{
  \code{\link{sw}}
  }
\examples{

## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## calculate robustness values
sw.rob(x)
  }
\keyword{misc}
