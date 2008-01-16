%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% on Tue Jan 15 18:33:06 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{normalizeWithinSlide.MAData}
\alias{normalizeWithinSlide.MAData}
\alias{MAData.normalizeWithinSlide}
\alias{normalizeWithinSlide.MAData}
\alias{normalizeWithinSlide,MAData-method}


 \title{Within-slide normalization}

 \description{
   Performs a within-slide normalization slide by slide.
   For a detailed explanation of normalization, see [1].\cr
 
   \emph{Note that
   the data in the object is replaced with the new normalized data and
   the old data is removed}. To keep the old data, make a copy of the
   object before normalizing by using \code{clone(ma)}, see
   \code{\link[R.oo]{clone.Object}} and example below.\cr
 }

 \usage{\method{normalizeWithinSlide}{MAData}(this, method, slides=NULL, weights=NULL, lowess=NULL, f=0.3, ...)}

 \arguments{
  \item{method}{The normalization method to be used. Currently there
   are four different methods;
   \code{"m"} - median normalization, which sets the median of log
   intensity ratios to zero,
   \code{"l"} - global lowess normalization,
   \code{"p"} - print-tip group lowess normalization, and
   \code{"s"} - scaled print-tip group lowess normalization.}
  \item{weights}{Weights between zero and one, that is, in [0,1], of each
   data point specifying how much that data points will affect the
   normalization. A data point with weight zero (or \code{\link[base:logical]{FALSE}}) will not affect
   the normalization, \emph{but} will be normalized. Currently only 0-1
   (or \code{\link[base:logical]{FALSE}}-\code{\link[base:logical]{TRUE}}) weights are supported. Non-zero weights are treated
   as ones.}
  \item{lowess}{When doing global lowess normalization, \code{method="l"},
   it is possible to specify the lowess line to be used. It is possible to
   specify individual lowess lines for each of the slides by letting
   lowess be a \emph{list} of lines. If \code{lowess=NULL} the lowess
   curve will be estimated from the data.}

   Note that \emph{only} one normalization is needed, i.e. doing different
   normalizations in serie on the same data set will not affect the
   results.

   Also note that it is only the log ratios, M, are affected by the
   normalization, i.e. the log intensities, A, are \emph{not} changed.
 }

 \references{
   \item{[1]}{S. Dudoit, Y. H. Yang, M. J. Callow, and T. P. Speed. Statistical
    methods for identifying differentially expressed genes in
    replicated cDNA microarray experiments (Statistics, UC Berkeley,
    Tech Report 578).
    URL: \url{http://www.stat.berkeley.edu/users/terry/zarray/Html/papersindex.html}}
 }

 \examples{
   # The option 'dataset' is used to annotate plots.
   options(dataset="sma:MouseArray")

   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   ma <- getSignal(raw)

   # Clone the data to get one non-normalized and one normalized data set.
   ma.norm <- clone(ma)

   # Normalize the data within slides using scaled print-tip normalization.
   normalizeWithinSlide(ma.norm, "s")

   # Plot data before and after normalization.
   layout(matrix(1:4, ncol=2, byrow=TRUE))
   plot(ma)
   plotSpatial(ma)
   plot(ma.norm)
   plotSpatial(ma.norm)
 }

 \note{
   Note that the layout must be set for print-tip (\code{method="p"}) and
   scaled (\code{method="s"}) normalization. If layout is not set, an
   exception will be thrown. Normally, the layout is already set, such as
   when the data is read from for instance GenePix, ScanAlyze and Spot.
 }

 \seealso{
   For across-slide normalization see \code{\link[aroma:normalizeAcrossSlides.MAData]{*normalizeAcrossSlides}()}.
   For more information see \code{\link{MAData}}.
 }

 \author{
   Henrik Bengtsson (\url{http://www.braju.com/R/}).
   Initial code for support of 'weights' by
   Jon McAuliffe, Statistics Dept, UC Berkeley.
   The original code was written by the sma authors
    Yee Hwa Yang \email{yeehwa@stat.berkeley.edu}
    Sandrine Dudoit \email{sandrine@stat.berkeley.edu} and
    Natalie Roberts \email{nroberts@wehi.edu.au}.
 }
\keyword{internal}
\keyword{methods}
