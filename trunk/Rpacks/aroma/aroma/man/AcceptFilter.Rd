%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AcceptFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AcceptFilter}
\docType{class}
\alias{AcceptFilter}


 \title{The AcceptFilter class}

 \description{
  Package:  aroma \cr
\bold{Class AcceptFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{AcceptFilter}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{GeneAcceptFilter}\cr

public static class \bold{AcceptFilter}\cr
extends \emph{\link[aroma]{SerialFilter}}\cr



   An AcceptFilter asks its input for indices and let only those indices
   through that is specified by its field \code{accept}.
 }

 \usage{AcceptFilter(input, accept=NULL, ...)}

 \arguments{
   \item{input}{The input \code{\link{MicroarrayData}} object.}
   \item{accept}{Vector of spot indices to be accepted.}
   \item{...}{Other arguments accepted by the constructor of the
    \code{\link{SerialFilter}} class.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:getIndex.AcceptFilter]{getIndex}} \tab Gets indices accepted by this filter.\cr
}


 \bold{Methods inherited from SerialFilter}:\cr
changeInput, getInput

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

   # Look at the top 5\% extreme M values
   fM  <- MFilter(tma, top=0.05, col="red")

   # And among those only look at the first 1000 spots.
   myFilter <- AcceptFilter(fM, accept=1:1000)

   plot(tma)
   highlight(myFilter, recursive=TRUE)
 }

 \seealso{
   See also the \code{\link{ParallelFilter}} class.
 }
\keyword{classes}
\keyword{internal}
\keyword{methods}
