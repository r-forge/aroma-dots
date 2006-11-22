%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ProbeLevelModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{fit.ProbeLevelModel}
\alias{fit.ProbeLevelModel}
\alias{ProbeLevelModel.fit}
\alias{fit.ProbeLevelModel}
\alias{fit,ProbeLevelModel-method}


 \title{Fits a unit groups (probeset) model}

 \description{
  Fits a unit groups (probeset) model to all or to a subset of the units.
  All estimates are stored to file.
  The non-array specific parameter estimates together with standard deviation
  estimates and convergence information are stored in one file.
  The parameter estimates specific to each array, typically "chip effects",
  are stored in array specific files.
 }

 \usage{\method{fit}{ProbeLevelModel}(this, units="remaining", ..., unitsPerChunk=moreUnits * 1e+05/length(getDataSet(this)), moreUnits=1, force=FALSE, verbose=FALSE)}

 \arguments{
   \item{...}{Not used.}
 }

 \value{
  Returns the indices of the units fitted, or \code{\link[base]{NULL}} if no units had to be fitted.
 }

 \details{
   Dataset-specific estimates [L = nbr of probes]:
    phi [L doubles] (probe affinities), sd(phi) [L doubles],
    isOutlier(phi) [L logicals]

   Algorithm-specific results:
    iter [1 integer], convergence1 [1 logical], convergence2 [1 logical]
    dTheta [1 double]
    sd(eps) - [1 double] estimated standard deviation of the error term

   Array-specific estimates [K = nbr of arrays]:
    theta [K doubles] (chip effects), sd(theta) [K doubles],
    isOutlier(theta) [K logicals]
 
   For each array and each unit group, we store:
     1 theta, 1 sd(theta), 1 isOutlier(theta), i.e. (float, float, bit)
   => For each array and each unit (with \eqn{G_j} groups), we store:
     \eqn{G_j} theta, \eqn{G_j} sd(theta), \eqn{G_j} isOutlier(theta),
   i.e. \eqn{G_j}*(float, float, bit).
   For optimal access we store all thetas first, then all sd(theta) and the
   all isOutlier(theta).
   To keep track of the number of groups in each unit, we have to have a
   (unit, ngroups) map.  This can be obtained from getUnitNames() for the
   AffymetrixCdfFile class.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{ProbeLevelModel}}.
 }


\keyword{internal}
\keyword{methods}
\keyword{IO}
