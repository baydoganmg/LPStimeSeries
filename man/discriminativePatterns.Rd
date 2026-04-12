\name{discriminativePatterns}
\alias{discriminativePatterns}
\title{Find and Visualize the Most Class-Discriminative Patterns}
\description{
  Identifies which terminal nodes (learned autopatterns) from a
  \code{learnPattern} ensemble are most discriminative for distinguishing
  between classes. Patterns are ranked by a chi-squared-like statistic
  that measures how much the per-class mean frequency in each terminal
  node deviates from the overall mean. The top patterns can be plotted
  with series colored by class label.
}
\usage{
discriminativePatterns(object, x, classes, n=5, plot=TRUE,
   orient=c(2,2), palette=NULL)
}
\arguments{
  \item{object}{an object of class \code{learnPattern}.}
  \item{x}{a matrix of time series (rows are series, columns are
    observations).}
  \item{classes}{a vector or factor of class labels, one per row of
    \code{x}.}
  \item{n}{number of top discriminative patterns to return (default 5).}
  \item{plot}{if \code{TRUE} (default), plot the top patterns colored by
    class.}
  \item{orient}{plot grid layout as \code{c(nrow, ncol)} (default
    \code{c(2,2)}).}
  \item{palette}{optional color palette vector. If \code{NULL}, colors
    are chosen automatically using \code{brewer.pal} or \code{rainbow}.}
}
\value{
  A list (returned invisibly) with the following components:
  \item{scores}{numeric vector of discriminability scores for all
    terminal nodes.}
  \item{ranking}{integer vector of terminal node indices ordered by
    decreasing discriminability.}
  \item{top}{integer vector of the top-\code{n} terminal node indices.}
  \item{class.means}{matrix of per-class mean frequencies (rows are
    classes, columns are terminal nodes).}
  \item{overall.means}{numeric vector of overall mean frequencies per
    terminal node.}
  \item{classes}{character vector of class level names.}
}
\note{
  The discriminability score for terminal node \eqn{j} is computed as
  \deqn{\sum_c n_c (\bar{f}_{cj} - \bar{f}_j)^2 / \bar{f}_j}
  where \eqn{n_c} is the number of series in class \eqn{c},
  \eqn{\bar{f}_{cj}} is the mean frequency of class \eqn{c} in node
  \eqn{j}, and \eqn{\bar{f}_j} is the overall mean frequency. Nodes
  where one class lands much more (or less) often than others receive
  high scores.

  When \code{plot=TRUE}, the predictor (triangle) and target (circle)
  segments are overlaid for all series, colored by class. At most
  \code{prod(orient)} patterns are plotted per page.
}
\references{
  Baydogan, M. G. and Runger, G. (2016), ``Time Series Representation and
  Similarity Based on Local Autopatterns'', Data Mining and Knowledge
  Discovery, 30(2), 476-509. \doi{10.1007/s10618-015-0425-y}.
}
\seealso{
  \code{\link{learnPattern}}, \code{\link{visualizePattern}},
  \code{\link{predict.learnPattern}}, \code{\link{plotMDS}}
}
\examples{
data(GunPoint)
set.seed(71)
ensemble <- learnPattern(GunPoint$trainseries)

## Find the 4 most discriminative patterns on the training set
dp <- discriminativePatterns(ensemble, GunPoint$trainseries,
                             GunPoint$trainclass, n=4)

## Inspect the top pattern scores
dp$scores[dp$top]
}
\author{Mustafa Gokce Baydogan}
\keyword{similarity}
\keyword{hplot}
