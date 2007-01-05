%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.extractSnpFeatureSet.R
% 
% on Sat Jan 06 01:47:04 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractSnpFeatureSet.AffymetrixCelSet}
\alias{extractSnpFeatureSet.AffymetrixCelSet}
\alias{AffymetrixCelSet.extractSnpFeatureSet}
\alias{extractSnpFeatureSet.AffymetrixCelSet}
\alias{extractSnpFeatureSet,AffymetrixCelSet-method}


\title{Extracts an in-memory SnpFeatureSet object from the CEL set}

\description{
 Extracts an in-memory SnpFeatureSet object from the CEL set.
 Note that any modifications done to the extract object will \emph{not}
 be reflected in the original CEL set.
}

\usage{\method{extractSnpFeatureSet}{AffymetrixCelSet}(this, ..., verbose=FALSE)}

\arguments{
  \item{...}{Argument passed to \code{oligo::read.cellfiles()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an \code{\link[oligo]{SnpFeatureSet-class}} object.
}

\details{
 Since the \pkg{oligo} package is making use of special CDF environment
 packages, this method will warn if the needed package is missing and
 explain that \pkg{oligo} will later try to download and install it
 automatically.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixCelSet}}.
}



\keyword{internal}
\keyword{methods}
\keyword{IO}
\keyword{programming}
