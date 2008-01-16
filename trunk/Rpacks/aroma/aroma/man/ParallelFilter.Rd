%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ParallelFilter.R
% 
% on Tue Jan 15 18:34:41 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{ParallelFilter}
\docType{class}
\alias{ParallelFilter}


\title{The ParallelFilter class}

\description{
 Package:  aroma \cr
\bold{Class ParallelFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{ParallelFilter}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{AndFilter}, \link[aroma]{OrFilter}\cr

public abstract static class \bold{ParallelFilter}\cr
extends \emph{\link[aroma]{Filter}}\cr



  A ParallelFilter is a filter that passes through indices from several
  inputs given some criteria.
  Examples of parallel filters are the \link{AndFilter} and the
  \code{\link{OrFilter}}, which provides the logical operators AND and OR,
  respectively.
}

\usage{ParallelFilter(..., cex=NULL, col=NULL, pch=NULL, visible=TRUE)}

\arguments{
  \item{...}{The input \code{\link{Filter}}s to be connected to.}
  \item{cex}{The scale factor of symbols that this filter highlights.}
  \item{col}{The color of symbols that this filter highlights.}
  \item{pch}{The plot symbols that this filter highlights.}
  \item{visible}{If \code{\link[base:logical]{TRUE}}, the data points filtered out by this filter
    will be highlighted, otherwise not.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab filters \tab The input filters. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:changeInput.ParallelFilter]{changeInput}} \tab Replaces input filters with other ones.\cr
 \tab \code{\link[aroma:getInput.ParallelFilter]{getInput}} \tab Gets the input filters.\cr
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

   plot(tma, "TvsSE");
   highlight(myFilter, recursive=TRUE);
}

\seealso{
  See also the \code{\link{SerialFilter}} class.
  For logical filters see \code{\link{AndFilter}}, \code{\link{OrFilter}}, and
  \code{\link{NotFilter}}.
  For data filters see \code{\link{MFilter}}, \code{\link{AFilter}}, \code{\link{TFilter}} and
  \code{\link{SEFilter}}.
}

\keyword{manip}
\keyword{classes}
