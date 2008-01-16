%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% on Tue Jan 15 18:35:58 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getBackground.RawData}
\alias{getBackground.RawData}
\alias{RawData.getBackground}
\alias{getBackground.RawData}
\alias{getBackground,RawData-method}


 \title{Gets the background signal}

 \description{
  Gets the background signal as a MAData object.
 }

 \usage{\method{getBackground}{RawData}(this, slides=NULL, ...)}

 \arguments{
   \item{slides}{Subset of slides to be returned. If \code{\link[base]{NULL}}, all slides
   are returned.}
 }

 \value{
   Returns the background signal as a \code{\link{MAData}} object.
 }

 \examples{
   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   bg <- getBackground(raw)
   plot(bg)
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma:getForeground.RawData]{*getForeground}()},
   \code{\link[aroma:getSignal.RawData]{*getSignal}()}.
   For more information see \code{\link{RawData}}.
 }
\keyword{internal}
\keyword{methods}
