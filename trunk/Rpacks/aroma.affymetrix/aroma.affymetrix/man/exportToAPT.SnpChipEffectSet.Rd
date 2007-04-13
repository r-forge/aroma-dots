%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.exportToAPT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{exportToAPT.SnpChipEffectSet}
\alias{exportToAPT.SnpChipEffectSet}
\alias{SnpChipEffectSet.exportToAPT}
\alias{exportToAPT.SnpChipEffectSet}
\alias{exportToAPT,SnpChipEffectSet-method}


\title{Export probeset summaries to an Affymetrix Power Tool (APT) result file}

\description{
 Export probeset summaries to an Affymetrix Power Tool (APT) result file named *.summary.txt.
}

\usage{\method{exportToAPT}{SnpChipEffectSet}(this, filename=sprintf("\%s,\%s.summary.txt", getFullName(this), getChipType(getCdf(this), fullname = FALSE)), path=NULL, sampleNames=NULL, ..., verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the APT result file.}
  \item{path}{An optional path to the file.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns (invisible) the pathname of the generated file.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SnpChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
