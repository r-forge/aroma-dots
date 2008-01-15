%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{convertFieldNames.GenePixData}
\alias{convertFieldNames.GenePixData}
\alias{GenePixData.convertFieldNames}
\alias{convertFieldNames.GenePixData}
\alias{convertFieldNames,GenePixData-method}


\title{Convert field names that contains patterns}

\description{
 Convert field names that contains patterns.
}

\usage{\method{convertFieldNames}{GenePixData}(this, patterns, channels=NULL, channelNames=NULL, ...)}

\arguments{
  \item{patterns}{A \code{\link[base]{vector}} of \code{\link[base]{character}} strings containing field names
    with patterns to be converted.}
  \item{channels}{An \code{\link[base]{integer}} \code{\link[base]{vector}}.}
  \item{channelNames}{A \code{\link[base]{character}} \code{\link[base]{vector}}.}
  \item{...}{Not used.}
}

\value{
 Returns the \code{\link[base]{vector}} of \code{\link[base]{character}} string of the converted field names.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  print(GenePixData$convertFieldNames(c("area", "F<channel[1]> Mean")));
}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
