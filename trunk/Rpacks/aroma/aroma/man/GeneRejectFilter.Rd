%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneRejectFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GeneRejectFilter}
\docType{class}
\alias{GeneRejectFilter}


\title{The GeneRejectFilter class}

\description{
 Package:  aroma \cr
\bold{Class GeneRejectFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma]{RejectFilter}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{GeneRejectFilter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{GeneRejectFilter}\cr
extends \link[aroma]{RejectFilter}\cr



  An GeneRejectFilter asks its input for indices and let only those
  indices through that do \emph{not} corresponds to a given set of genes.
}

\usage{GeneRejectFilter(input, layout, genes=NULL, ids=NULL, names=NULL, ...)}

\arguments{
  \item{input}{The input \code{\link{Filter}} to be connected to.}
  \item{layout}{A \code{\link{Layout}} object.}
  \item{genes}{Gene}.
  \item{ids}{}.
  \item{names}{}.
  \item{...}{Any arguments accepted by the \code{\link{RejectFilter}} constructor.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getIndex} \tab  -\cr
}


 \bold{Methods inherited from RejectFilter}:\cr
getIndex

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

   # However, the blanks, e.g. empty spots, are not of interest
   myFilter <- GeneRejectFilter(fM, layout=layout, ids="BLANK", col="orange")

   plot(tma)
   highlight(myFilter, recursive=TRUE)
}

\seealso{
  See also the \code{\link{RejectFilter}} class.
}

\keyword{manip}
\keyword{classes}
