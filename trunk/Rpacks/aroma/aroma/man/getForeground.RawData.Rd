%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RawData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getForeground.RawData}
\alias{getForeground.RawData}
\alias{RawData.getForeground}
\alias{getForeground.RawData}
\alias{getForeground,RawData-method}


 \title{Gets the foreground signal}

 \description{
  Gets the foreground signal as a MAData object.
 }

 \usage{\method{getForeground}{RawData}(this, slides=NULL, ...)}

 \arguments{
   \item{slides}{Subset of slides to be returned. If \code{\link[base]{NULL}}, all slides
   are returned.}
 }

 \value{
   Returns the foreground signal as a \code{\link{MAData}} object.
 }

 \examples{
   SMA$loadData("mouse.data")
   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
   raw <- RawData(mouse.data, layout=layout)
   fg <- getForeground(raw)
   plot(fg)
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma:getBackground.RawData]{*getBackground}()},
   \code{\link[aroma:getSignal.RawData]{*getSignal}()}.
   For more information see \code{\link{RawData}}.
 }
\keyword{internal}
\keyword{methods}
