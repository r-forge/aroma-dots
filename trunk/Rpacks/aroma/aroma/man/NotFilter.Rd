%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  NotFilter.R
% 
% on Tue Jan 15 18:34:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{NotFilter}
\docType{class}
\alias{NotFilter}


\title{The NotFilter class}

\description{
 Package:  aroma \cr
\bold{Class NotFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{NotFilter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{NotFilter}\cr
extends \emph{\link[aroma]{SerialFilter}}\cr



  A NotFilter takes a set of indices from a single input and returns
  the complementary set of indices.
}

\usage{NotFilter(input, ...)}

\arguments{
  \item{input}{The input \code{\link{Filter}} to be connected to.}
  \item{...}{Any arguments accepted by the \code{\link{Filter}} constructor.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getIndex} \tab  -\cr
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

   fM  <- MFilter(tma, top=0.05, col="red")
   fT  <- TFilter(tma, top=0.05, col="blue")
   fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
   fSE <- NotFilter(fNotSE, visible=FALSE)
   myFilter <- AndFilter(fM, fT, fSE, col="purple")

   plot(tma, "TvsSE");
   highlight(myFilter, recursive=TRUE);
}

\seealso{
  See also the \code{\link{ParallelFilter}} class.
}

\keyword{manip}
\keyword{classes}
