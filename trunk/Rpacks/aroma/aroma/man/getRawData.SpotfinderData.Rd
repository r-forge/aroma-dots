%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotfinderData.R
% 
% on Tue Jan 15 18:37:29 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.SpotfinderData}
\alias{getRawData.SpotfinderData}
\alias{SpotfinderData.getRawData}
\alias{getRawData.SpotfinderData}
\alias{getRawData,SpotfinderData-method}


\title{Gets the raw intensites from the SpotfinderData structure}

\description{
 Extracts the red and green spot intensitites (both foreground and background)
 from the SpotfinderData object and returns a \code{\link{RawData}} object.
}

\usage{\method{getRawData}{SpotfinderData}(this, slides=NULL, ...)}

\arguments{
  \item{slides}{Specifying which slides to be extracted. If \code{\link[base]{NULL}},
    all slides are considered.}
}

\value{
  Returns a \code{\link{RawData}} object containing the specified slides.
  Note that the returned foreground and background signals for
  Spotfinder data is \emph{not} necessarily \code{\link[base]{integer}}s,
  e.g. \code{R==0.32} is possible.
}

\details{
  The R and Rb channels will come from the *A* fields, and
  the G and Gb channels will come from the *B* fields.
  To swap the channels use dyeSwap().
}

\examples{\dontrun{# For an example see help(SpotfinderData).}}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SpotfinderData}}.
}
\keyword{internal}
\keyword{methods}
