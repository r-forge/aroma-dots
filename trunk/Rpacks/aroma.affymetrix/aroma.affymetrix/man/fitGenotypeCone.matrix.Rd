%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fitGenotypeCone.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{fitGenotypeCone.matrix}
\alias{fitGenotypeCone.matrix}
\alias{fitGenotypeCone.matrix}


 \title{Fits an affine transformation to allele A and allele B data}

 \description{
  Fits an affine transformation to allele A and allele B data using robust estimators.
 }

 \usage{\method{fitGenotypeCone}{matrix}(y, alpha=c(0.1, 0.075, 0.05, 0.03, 0.01), q=2, Q=98, ...)}

 \arguments{
   \item{y}{A \code{\link[base]{numeric}} Nx2 \code{\link[base]{matrix}} with one column for each allele and
      where N is the number of data points.}
   \item{alpha}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of decreasing values in (0,1).
      This parameter "determines how far we are willing to press the
      boundary of the [genotype cone]".  Lowering \code{alpha} expand
      the cone.  When \code{alpha} goes to zero, all data points will
      be on or inside the cone.}
   \item{q,Q}{Percentiles in [0,100] for which data points that are
      below (above) will be assigned zero weight in the fitting of
      the parameters.}
   \item{...}{Additional arguments passed to \code{\link[sfit]{cfit}}.}
 }

 \value{
   Returns the parameter estimates as a named \code{\link[base]{list}} with elements:
   \itemize{
    \item{M}{An estimate of the three vertices defining the genotype
      triangle.  These three vertices are describes as an 2x3 \code{\link[base]{matrix}}
      with column \code{origin}, \code{AA}, and \code{BB}.}
    \item{Minv}{The inverse of \code{M}.}
    \item{origin}{The estimate of the shift.}
    \item{W}{The estimate of shear/rotation matrix with columns
             \code{AA} and \code{BB}.}
    \item{Winv}{The inverse of \code{W}.}
    \item{params}{The parameters used for the fit, i.e.
       \code{alpha}, \code{q}, \code{Q}, and  those passed in \code{...}.}
   }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
  To backtransform data fitted using this method,
  see \code{\link{backtransformGenotypeCone}}().
  Internally \code{\link[sfit]{cfit}} of the \pkg{sfit} package is used.
 }


\keyword{methods}
\keyword{internal}
