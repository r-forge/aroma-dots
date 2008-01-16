%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% on Tue Jan 15 18:35:58 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getColors.RawData}
\alias{getColors.RawData}
\alias{RawData.getColors}
\alias{getColors.RawData}
\alias{getColors,RawData-method}


 \title{Generates colors for each of the specified spots}

 \description{
  Generates colors for each of the specified spots, which can be passed to
  the \code{col} argument in most plot functions. However, note that most
  of the plot functions in this package automatically recognize palettes
  if they are passed to the \code{col} argument.
 }

 \usage{\method{getColors}{RawData}(this, what=c("R", "G"), slide=1, include=NULL, exclude=NULL, palette=NULL, bgSubtract=FALSE, bg.subtract=bgSubtract, range=matrix(c(0, 16, 0, 16), nrow = 2), ...)}

 \arguments{
  \item{what}{The fields to be used.}
  \item{slide}{Specifies for which slide the colors should be generated for.}
  \item{include}{The indices of the spots that should be included.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are considered.
   If \code{\link[base]{NULL}} all spots are considered.}
  \item{exclude}{The indices of the spots that should be excluded.
   If it is instead a name of one or more flags, the spots which have been
   flagged with these flags are excluded.}
   \item{bgSubstract}{If \code{\link[base:logical]{TRUE}}, the background is subtracted from the
         foreground, before the colors are generated.}
   \item{palette}{The palette to be used. If \code{\link[base]{NULL}} or \code{"auto"}
         the palette that is most applicable to the value of \code{what}
         will be used.}
   \item{range}{Cut-off range for data values. The format is a matrix with
     two rows and \code{length(what)} number of columns.}
 }

 \value{Returns a \code{\link[base]{vector}} of colors.}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \examples{
   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   col <- getColors(raw, what="G", slide=4);
 }

 \seealso{
   For more information see \code{\link{RawData}}.
 }
\keyword{internal}
\keyword{methods}
