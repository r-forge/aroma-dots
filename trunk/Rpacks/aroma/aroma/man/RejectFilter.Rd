%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RejectFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{RejectFilter}
\docType{class}
\alias{RejectFilter}


\title{The RejectFilter class}

\description{
 Package:  aroma \cr
\bold{Class RejectFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{RejectFilter}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{GeneRejectFilter}\cr

public static class \bold{RejectFilter}\cr
extends \emph{\link[aroma]{SerialFilter}}\cr



  An RejectFilter asks its input for indices and let only those indices
  through that is specified by its field \code{reject}.
}

\usage{RejectFilter(input, reject=NULL, ...)}

\arguments{
  \item{input}{The input \code{\link{Filter}} to be connected to.}
  \item{reject}{The indices to be rejected by this filter.}
  \item{...}{Any arguments accepted by the \code{\link{SerialFilter}} constructor.}
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
$, $<-, [[, [[<-, as.character, attach, attachLocally, clone, detach, environment, equals, extend, finalize, getFields, getInstanciationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save

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

   # However, the first 1000 spots are not of interest
   myFilter <- RejectFilter(fM, reject=1:1000)

   plot(tma);
   highlight(myFilter, recursive=TRUE);
}

\seealso{
  See also the \code{\link{ParallelFilter}} class.
}

\keyword{manip}
\keyword{classes}
