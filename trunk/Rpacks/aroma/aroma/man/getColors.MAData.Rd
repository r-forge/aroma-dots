%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getColors.MAData}
\alias{getColors.MAData}
\alias{MAData.getColors}
\alias{getColors.MAData}
\alias{getColors,MAData-method}


 \title{Generates colors for each of the specified spots}

 \description{
  Generates colors for each of the specified spots, which can be passed as
  a \code{col} argument in most plot functions.
 }

 \usage{\method{getColors}{MAData}(this, what=c("A", "M"), slide=1, include=NULL, exclude=NULL, palette=NULL, range="auto", ...)}

 \arguments{
  \item{slide}{Specifies for which slide the colors should be generated for.}
  \item{include}{The spot indices to be plotted. If \code{\link[base]{NULL}} all spots are considered.}
  \item{exclude}{The spot indices \emph{not} to be plotted. If \code{\link[base]{NULL}} no spots are excluded (from all or from the one specified by \code{include}).}
  \item{palette}{The palette to be used. Currently, only the \code{"redgreen"} palette is available.}
  \item{M.range}{The range of "normal" M values. All M values outside this range will be cutoff.}
  \item{A.range}{The range of "normal" A values. All A values outside this range will be cutoff.}
 }

 \value{Returns a \code{\link[base]{vector}} of colors.}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \examples{
   # The option 'dataset' is used to annotate plots.
   options(dataset="sma:MouseArray")

   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   ma <- getSignal(raw, bgSubtract=TRUE)
   col <- getColors(ma, slide=4)
 }

 \seealso{
   \code{\link[R.graphics]{Colors}}.
   For more information see \code{\link{MAData}}.
 }
\keyword{internal}
\keyword{methods}
