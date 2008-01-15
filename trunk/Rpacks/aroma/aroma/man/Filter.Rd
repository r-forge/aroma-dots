%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Filter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Filter}
\docType{class}
\alias{Filter}


\title{The Filter class}

\description{
 Package:  aroma \cr
\bold{Class Filter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{Filter}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{AcceptFilter}, \link[aroma]{AFilter}, \link[aroma]{AndFilter}, \link[aroma]{BFilter}, \link[aroma]{DfFilter}, \link[aroma]{FieldFilter}, \link[aroma]{GeneAcceptFilter}, \link[aroma]{GeneRejectFilter}, \link[aroma]{MFilter}, \link[aroma]{NotFilter}, \link[aroma]{OrFilter}, \emph{\link[aroma]{ParallelFilter}}, \link[aroma]{RejectFilter}, \link[aroma]{SEFilter}, \emph{\link[aroma]{SerialFilter}}, \link[aroma]{TFilter}\cr

public abstract static class \bold{Filter}\cr
extends \link[R.oo]{Object}\cr



  A Filter object takes a vector of indices from its input(s) and passes
  some or all of them through given some criteria. There are two main
  classes of filters; (1) the \code{\link{SerialFilter}} class and (2) the
  \code{\link{ParallelFilter}} class. Shortly, the \emph{serial} filters gets
  some indices from one single source, whereas the \emph{parallel}
  filters gets indices from several sources. Example of serial filters
  are the \code{\link{FieldFilter}}s and the \code{\link{NotFilter}}. Example of
  parallel filters are the \emph{logical} filters \code{\link{AndFilter}} and
  \code{\link{OrFilter}}.
}

\usage{Filter(cex=NULL, col=NULL, pch=20, visible=TRUE)}

\arguments{
  \item{cex}{The scale factor of symbols that this filter highlights.}
  \item{col}{The color of symbols that this filter highlights.}
  \item{pch}{The plot symbols that this filter highlights.}
  \item{visible}{If \code{\link[base:logical]{TRUE}}, the data points filtered out by this filter
    will be highlighted, otherwise not.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab parameter \tab A list of parameters for the filter. \cr
   \tab visible \tab Specifies if the data through this filter should be
        displayed in plots etc or not. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{\link[aroma:changeInput.Filter]{changeInput}} \tab Change input(s) on this filter and optionally all filters down the stream.\cr
 \tab \code{\link[aroma:getIndex.Filter]{getIndex}} \tab Gets the indices accepted by the filter.\cr
 \tab \code{\link[aroma:getInput.Filter]{getInput}} \tab Gets all the input objects connected to the filter.\cr
 \tab \code{\link[aroma:getParameter.Filter]{getParameter}} \tab Gets the values for a specific parameter for indices accepted by the filter.\cr
 \tab \code{\link[aroma:getVisible.Filter]{getVisible}} \tab Checks if the filter is visible or not.\cr
 \tab \code{\link[aroma:highlight.Filter]{highlight}} \tab Highlights the data points accepted by the filter.\cr
 \tab \code{\link[aroma:setParameter.Filter]{setParameter}} \tab Sets the values for a specific parameter for indices accepted by the filter.\cr
 \tab \code{\link[aroma:setVisible.Filter]{setVisible}} \tab Sets if the filter should bee visible or not.\cr
 \tab \code{\link[aroma:text.Filter]{text}} \tab Labels the data points accepted by the filter.\cr
}


\bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)

   ma <- getSignal(raw)
   normalizeWithinSlide(ma, method="s")
   normalizeAcrossSlides(ma)

   tma <- as.TMAData(ma)

   # Filter out too weak spots, i.e. A < 8 (<256).
   fA  <- AFilter(tma, range=c(8, Inf), visible=FALSE)

   # Look at the top 5\% M-values.
   fM  <- MFilter(tma, top=0.05, col="red")

   # Look at the top 5\% T-values, but not those with to small SE.
   fT  <- TFilter(tma, top=0.05, col="blue")
   fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
   fSE <- NotFilter(fNotSE, visible=FALSE)

   # Selects those spot that passes through all filters.
   myFilter <- AndFilter(fA, fM, fT, fSE, col="purple", pch=19)

   plot(tma, "TvsSE", xlog=2)
   highlight(myFilter, recursive=TRUE)
}

\keyword{manip}
\keyword{classes}
