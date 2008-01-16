%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SerialFilter.R
% 
% on Tue Jan 15 18:37:04 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SerialFilter}
\docType{class}
\alias{SerialFilter}


\title{The SerialFilter class}

\description{
 Package:  aroma \cr
\bold{Class SerialFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{SerialFilter}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{AcceptFilter}, \link[aroma]{AFilter}, \link[aroma]{BFilter}, \link[aroma]{DfFilter}, \link[aroma]{FieldFilter}, \link[aroma]{GeneAcceptFilter}, \link[aroma]{GeneRejectFilter}, \link[aroma]{MFilter}, \link[aroma]{NotFilter}, \link[aroma]{RejectFilter}, \link[aroma]{SEFilter}, \link[aroma]{TFilter}\cr

public abstract static class \bold{SerialFilter}\cr
extends \emph{\link[aroma]{Filter}}\cr



  A SerialFilter is a filter that passes through indices from a single
  input given some criteria.
}

\usage{SerialFilter(input=NULL, ...)}

\arguments{
  \item{input}{The input \code{\link{Filter}} to be connected to.}
  \item{...}{Any arguments accepted by the \code{\link{Filter}} constructor.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab input \tab The input. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:changeInput.SerialFilter]{changeInput}} \tab Change input(s) on this filter and optionally all filters down the stream.\cr
 \tab \code{\link[aroma:getInput.SerialFilter]{getInput}} \tab Gets all the input objects connected to the filter.\cr
}


 \bold{Methods inherited from Filter}:\cr
as.character, changeInput, getIndex, getInput, getParameter, getVisible, highlight, setParameter, setVisible, text

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save


}


\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)

   ma <- getSignal(raw)
   normalizeWithinSlide(ma, "s")
   normalizeAcrossSlides(ma)

   tma <- as.TMAData(ma)

   fM  <- MFilter(tma, top=0.05, col="red")
   fT  <- TFilter(tma, top=0.05, col="blue")
   fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
   fSE <- NotFilter(fNotSE, visible=FALSE)
   myFilter <- AndFilter(fM, fT, fSE, col="purple")

   plot(tma, "TvsSE")
   highlight(myFilter, recursive=TRUE)
}

\seealso{
  See also the \code{\link{ParallelFilter}} class.
}
\keyword{classes}
