%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.extractSnpQSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{extractSnpQSet.SnpChipEffectSet}
\alias{extractSnpQSet.SnpChipEffectSet}
\alias{SnpChipEffectSet.extractSnpQSet}
\alias{extractSnpQSet.SnpChipEffectSet}
\alias{extractSnpQSet,SnpChipEffectSet-method}


\title{Extracts chip effects as a SnpQSet}

\description{
 Extracts chip effects as a SnpQSet.
 The chip effects are extracted to the log (base 2) scale.
 Note: All chip effects are copied into memory.  Moreover, any modifications
 of the extracted data is \emph{not} reflected in the data files, and vice
 versa.
}

\usage{\method{extractSnpQSet}{SnpChipEffectSet}(this, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[oligo]{SnpQSet-class}} object.
}

\details{
  Simple benchmarking on shadowfax:
  To extract chip effects for the 90 CEPH Mapping50K\_Xba240 arrays,
  it takes ~90 seconds, i.e. 1s/array.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SnpChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
