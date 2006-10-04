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

\usage{\method{extractSnpQSet}{SnpChipEffectSet}(this, ..., verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns a \code{\link[oligo]{SnpQSet-class}} object.
}

\details{
  Simple benchmarking on shadowfax: To fit crlmm() on 90 CEPH 50K chips
  it takes ~100 minutes, i.e. 67 s/array.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SnpChipEffectSet}}.
}
\keyword{internal}
\keyword{methods}
