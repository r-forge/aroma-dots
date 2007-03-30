%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeFragmentLength.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{normalizeFragmentLength}
\alias{normalizeFragmentLength.default}
\alias{normalizeFragmentLength}

 \title{Normalizes signals for PCR fragment-length effects}

 \description{
  Normalizes signals for PCR fragment-length effects. Some or all signals are used to estimated the
  normalization function.  All signals are normalized.
 }

 \usage{\method{normalizeFragmentLength}{default}(y, fragmentLengths, targetFcn=NULL, subsetToFit=NULL, ...)}

 \arguments{
   \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length K of signals to be normalized}
   \item{fragmentLengths}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of length K of fragment lengths.}
   \item{targetFcn}{}
   \item{subsetToFit}{The subset of data points used to fit the
      normalization function.
      If \code{\link[base]{NULL}}, all data points are considered.}
   \item{...}{Additional arguments passed to \code{\link[stats]{lowess}}.}
 }

 \value{
   Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of the normalized signals.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
