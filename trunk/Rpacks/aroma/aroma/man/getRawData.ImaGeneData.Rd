%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ImaGeneData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.ImaGeneData}
\alias{getRawData.ImaGeneData}
\alias{ImaGeneData.getRawData}
\alias{getRawData.ImaGeneData}
\alias{getRawData,ImaGeneData-method}


\title{Gets the raw intensites from two ImaGeneData objects}

\description{
 Extracts the red and green spot intensitites (both foreground and
 background) from two ImaGeneData objects and returns a
 \code{\link{RawData}} object.
}

\usage{\method{getRawData}{ImaGeneData}(igR, igG, slides=NULL, fg="median", bg="mean", ...)}

\arguments{
  \item{igR, igG}{The "red" and the "green" ImageGeneData objects.}
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
  The R and Rb channels will come from the igR object, and
  the G and Gb channels will come from the igG object.
  To swap the channels swap the arguments or use dyeSwap() afterwards.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ImaGeneData}}.
}
\keyword{internal}
\keyword{methods}
