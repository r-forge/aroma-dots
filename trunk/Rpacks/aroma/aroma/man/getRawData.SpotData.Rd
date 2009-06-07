%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.SpotData}
\alias{getRawData.SpotData}
\alias{SpotData.getRawData}
\alias{getRawData.SpotData}
\alias{getRawData,SpotData-method}


\title{Gets the raw intensites from the SpotData structure}

\description{
 Extracts the red and green spot intensitites (both foreground and background)
 from the SpotData object and returns a \code{\link{RawData}} object.
}

\usage{\method{getRawData}{SpotData}(this, slides=NULL, fgs="auto", bgs="auto", ...)}

\arguments{
  \item{slides}{Specifying which slides to be extracted. If \code{\link[base]{NULL}},
    all slides are considered.}
  \item{fg}{If \code{"mean"}, the mean foreground intensities are returned.
    If \code{"median"}, the median foreground intensities are returned.}
  \item{bg}{If \code{"mean"}, the mean background intensities are returned.
    If \code{"median"}, the median background intensities are returned. To
    get the mean of the median of the background intensites in the four
    diamond shaped ares around the spots, use the value \code{"valley"}.
    There are also three morphological background estimates; \code{morph},
    \code{morph.erode}, and \code{morph.close.open}.}
}

\value{
  Returns a \code{\link{RawData}} object containing the specified slides.
}

\details{
  The R and Rb channels will come from the *R* fields, and
  the G and Gb channels will come from the *G* fields.
  To swap the channels use dyeSwap().
}

\examples{\dontrun{# For an example see help(SpotData).}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
