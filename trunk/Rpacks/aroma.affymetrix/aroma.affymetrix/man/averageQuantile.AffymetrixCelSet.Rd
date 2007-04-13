%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{averageQuantile.AffymetrixCelSet}
\alias{averageQuantile.AffymetrixCelSet}
\alias{AffymetrixCelSet.averageQuantile}
\alias{averageQuantile.AffymetrixCelSet}
\alias{averageQuantile,AffymetrixCelSet-method}


 \title{Gets the average empirical distribution across all samples}

 \description{
  Gets the average empirical distribution across all samples.
 }

 \usage{\method{averageQuantile}{AffymetrixCelSet}(this, probes=NULL, excludeCells=NULL, ..., verbose=FALSE)}

 \arguments{
   \item{probes}{An optional \code{\link[base]{numeric}} \code{\link[base]{vector}} specifying what subset of
      probes to be used to calculate the empirical distribution.
      If \code{\link[base]{NULL}}, all probes are used.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
  Returns a \code{\link[base]{double}} \code{\link[base]{vector}} of length N.
 }

 \section{Missing values}{
   If \code{\link[base]{NA}}s are detected in a sample, these are excluded and the
   \code{approx()} function (\code{\link[stats]{approx}}) is used to "expand"
   the \code{\link[base]{vector}} of the remaining values so that the sorted \code{\link[base]{vector}}
   (still) have length N.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.light]{averageQuantile.list}}
   For more information see \code{\link{AffymetrixCelSet}}.
 }
\keyword{internal}
\keyword{methods}
