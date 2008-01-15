%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizePlatewise.MicroarrayData}
\alias{normalizePlatewise.MicroarrayData}
\alias{MicroarrayData.normalizePlatewise}
\alias{normalizePlatewise.MicroarrayData}
\alias{normalizePlatewise,MicroarrayData-method}


 \title{Normalization performed plate by plate}

 \description{
   Performs a normalization plate by plate on each of the slide seperately.
   For details on plate-wise normalization see [1] and for details
   on intensity dependent normalization see [2].\cr
 
   \emph{Note that
   the data in the object is replaced with the new normalized data and
   the old data is removed}. To keep the old data, make a copy of the
   object before normalizing by using \code{clone(ma)}, see
   \code{\link[R.oo]{clone.Object}} and example below.\cr
 }

 \usage{\method{normalizePlatewise}{MicroarrayData}(this, field, method="constant", slides=NULL, ...)}

 \arguments{
  \item{field}{The data field to be normalized.}
  \item{method}{The normalization method to be used. If \code{"constant"}
   the values of the data field will be shifted to have zero bias within
   each plate group. If \code{"A"}, the values within each plate group
   will be normalized against intensity dependent effects.}
  \item{slides}{Slides to be included in the normalization. If \code{\link[base]{NULL}},
     all slides are normalized.}
  \item{...}{Other arguments accepted by underlying normalization methods.}
 }

 \details{
  The data will be normalize within each plate group \emph{individually}
  by either a) assuming same shift and scale for whole plate i) shifting
  the data to have median/mean zero. ii) rescale the data so all the
  groups have the same median absolute deviation (MAD) or standard
  deviation (sd) or by b) assuming intensity dependent effects. For
  more details see [1].
 }

 \examples{
   # The option 'dataset' is used to annotate plots.
   options(dataset="sma:MouseArray")

   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   ma <- getSignal(raw)

   # Scaled intensity normalization print-tip by print-tip
   ma.norm1 <- clone(ma)
   normalizeWithinSlide(ma.norm1, method="s")

   # Intensity normalization plate by plate
   ma.norm2 <- clone(ma.norm1)
   normalizePlatewise(ma.norm2, field="M", method="A")

   # Plot data before and after normalization.
   layout(matrix(1:9, ncol=3, byrow=TRUE))
   plot(ma)
   plotSpatial(ma)
   plotPrintorder(ma)
   plot(ma.norm1)
   plotSpatial(ma.norm1)
   plotPrintorder(ma.norm1)
   plot(ma.norm2)
   plotSpatial(ma.norm2)
   plotPrintorder(ma.norm2)
 }

 \references{
  \item{[1]}{Henrik Bengtsson, Plate Effects in cDNA microarray data,
    Matemathical Statistics, Centre for Matematical Sciences,
    Lund University, Sweden. Manuscript, 2002.}
  \item{[2]}{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
    methods for identifying differentially expressed genes in
    replicated cDNA microarray experiments (Statistics, UC Berkeley,
    Tech Report 578).
    URL: \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html}}
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For other within-slide normalization see \code{\link{MAData.normalizeWithinSlide}}.
   For across-slide normalization see \code{\link{MAData.normalizeAcrossSlides}}.
   For more information see \code{\link{MicroarrayData}}.
 }
\keyword{internal}
\keyword{methods}
