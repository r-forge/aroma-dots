%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.QuantArrayData}
\alias{getRawData.QuantArrayData}
\alias{QuantArrayData.getRawData}
\alias{getRawData.QuantArrayData}
\alias{getRawData,QuantArrayData-method}


\title{Gets the raw intensites from the QuantArray data structure}

\description{
 Extracts the red and green spot intensitites (both foreground and background)
 from the QuantArrayData object and returns a \code{\link{RawData}} object.
}

\usage{\method{getRawData}{QuantArrayData}(qa, slides=NULL, fg="mean", bg="mean", ...)}

\arguments{
  \item{slides}{Specifying which slides to be extracted. If \code{\link[base]{NULL}},
    all slides are considered.}
  \item{fg}{If \code{"mean"}, the mean foreground intensities are returned.
    If \code{"median"}, the median foreground intensities are returned.}
  \item{bg}{If \code{"mean"}, the mean background intensities are returned.
    If \code{"median"}, the median background intensities are returned.}
}

\value{
  Returns a \code{\link{RawData}} object containing the specified slides.
}

\details{
  The R and Rb channels will come from the ch1* fields, and
  the G and Gb channels will come from the ch2* fields.
  To swap the channels just use dyeSwap().
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{QuantArrayData}}.
}
\keyword{internal}
\keyword{methods}
