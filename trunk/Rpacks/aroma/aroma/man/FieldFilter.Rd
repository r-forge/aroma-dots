%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FieldFilter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FieldFilter}
\docType{class}
\alias{FieldFilter}


\title{The FieldFilter class}

\description{
 Package:  aroma \cr
\bold{Class FieldFilter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma]{Filter}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma]{SerialFilter}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FieldFilter}\cr

\bold{Directly known subclasses:}\cr
\link[aroma]{AFilter}, \link[aroma]{BFilter}, \link[aroma]{DfFilter}, \link[aroma]{MFilter}, \link[aroma]{SEFilter}, \link[aroma]{TFilter}\cr

public static class \bold{FieldFilter}\cr
extends \emph{\link[aroma]{SerialFilter}}\cr



  A FieldFilter provides methods to extract indices from
  \code{\link{MicroarrayData}} objects given some criteria on one of the
  \emph{fields}. Possible criterias on this field are \emph{top},
  and \emph{range}.
  The criteria \emph{top} filters out the top \emph{N} fraction
  (or number if \emph{N>1}).
  The criteria \emph{range} filters accepts the spots with a field value
  within the given range.
}

\usage{FieldFilter(mad, field, bottom=NULL, top=NULL, range=NULL, absolute.values=TRUE, cex=NULL, col=NULL, pch=NULL, visible=TRUE)}

\arguments{
  \item{mad}{A \code{\link{MicroarrayData}} object to be filtered.}
  \item{field}{The field (\code{\link[base]{character}} name) to be filtered.}
  \item{top, bottom}{If specified, to top (bottom) values are filtered out.
   If an \code{\link[base]{integer}} one or greater, that number of indicies will be passed.
   If a \code{\link[base]{numeric}} between zero and one, that ratio will be passed.}
  \item{range}{The range of values to be passed.}
  \item{absolute.values}{If \code{\link[base:logical]{TRUE}}, absolute values are filtered,
   otherwise not.}
  \item{cex}{The scale factor of symbols that this filter highlights.}
  \item{col}{The color of symbols that this filter highlights.}
  \item{pch}{The plot symbols that this filter highlights.}
  \item{visible}{If \code{\link[base:logical]{TRUE}}, the data points filtered out by this filter
    will be highlighted, otherwise not.}
}

\section{Fields and Methods}{
 \bold{Fields}
 \tabular{rll}{
   \tab range \tab The range criteria. \cr
   \tab top \tab The top criteria. \cr
 }

 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma:getIndex.FieldFilter]{getIndex}} \tab Gets indices accepted by this filter.\cr
 \tab \code{getIndex0} \tab  -\cr
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

   fM  <- MFilter(tma, top=0.05, col="red")
   fT  <- TFilter(tma, top=0.05, col="blue")
   fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
   fSE <- NotFilter(fNotSE, visible=FALSE)
   myFilter <- AndFilter(fM, fT, fSE, col="purple")

   plot(tma, "TvsSE");
   highlight(myFilter, recursive=TRUE);
}

\seealso{
  See also the \code{\link{Filter}} class.
  For logical filters see \code{\link{AndFilter}}, \code{\link{OrFilter}}, and
  \code{\link{NotFilter}}.
  For data filters see \code{\link{MFilter}}, \code{\link{AFilter}}, \code{\link{TFilter}} and
  \code{\link{SEFilter}}.
}

\keyword{manip}
\keyword{classes}
