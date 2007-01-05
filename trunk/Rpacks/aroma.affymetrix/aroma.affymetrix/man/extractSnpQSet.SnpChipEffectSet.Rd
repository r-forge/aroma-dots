%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SnpChipEffectSet.extractSnpQSet.R
% 
% on Sat Jan 06 01:57:23 2007.
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

\usage{\method{extractSnpQSet}{SnpChipEffectSet}(this, units=NULL, transform=c("log", "asinh"), naValue=-5, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{units}{Units to be extracted.  If \code{\link[base]{NULL}}, all units are considered.
   Non-SNP units are to extracted.}
  \item{transform}{A \code{\link[base]{character}} string specifying what type of
   transformation should be applied to the chip-effect estimates, which
   are on the intensity scale.}
  \item{naValue}{A \code{\link[base]{double}} value used to replace \code{\link[base]{NA}} values.}
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached values are ignored.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
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
