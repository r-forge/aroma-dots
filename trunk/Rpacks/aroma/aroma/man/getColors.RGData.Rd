%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.R
% 
% on Tue Jan 15 18:36:16 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getColors.RGData}
\alias{getColors.RGData}
\alias{RGData.getColors}
\alias{getColors.RGData}
\alias{getColors,RGData-method}


 \title{Generates red to green colors for each of the specified spots}

 \description{
   Generates red to green colors for each of the specified spots, which can
  be passed to the \code{col} argument in most plot functions.
 }

 \usage{\method{getColors}{RGData}(this, what=c("A", "M"), slide=1, include=NULL, exclude=NULL, palette=NULL, range=matrix(c(0, 16, -2, 2), nrow = 2), ...)}

 \arguments{
  \item{palette}{The color palette to be used.}
  \item{M.range}{The range of the M=log2(R/G). Values outside this range will be saturated to the extreme green and red, respectively.}
  \item{A.range}{The range of the A=1/2*log2(R*G). Values outside this range will be saturated to the dark and light, respectively.}
 }

 \value{Returns nothing.}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \examples{
   SMA$loadData(c("mouse.data", "mouse.setup"))
   raw <- RawData(mouse.data, layout=as.Layout(mouse.setup))
   ma <- getSignal(raw)
   rg <- as.RGData(ma)
   col <- getColors(rg)
 }

 \seealso{
   For more information see \code{\link{RGData}}.
 }
\keyword{internal}
\keyword{methods}
