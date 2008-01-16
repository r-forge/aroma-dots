%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  RGData.fitIWPCA.R
% 
% on Tue Jan 15 18:36:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{fitMultiscanAffine.RGData}
\alias{fitMultiscanAffine.RGData}
\alias{RGData.fitMultiscanAffine}
\alias{fitMultiscanAffine.RGData}
\alias{fitMultiscanAffine,RGData-method}


 \encoding{latin1}

 \title{Fits an affine model to signals from a multi-scanned slide}

 \description{
   Fits an affine model to signals from a multi-scanned slide.

   This function is used internally by \code{\link[aroma:calibrateMultiscan.RGData]{*calibrateMultiscan}()} and
   is not intended to be used by the end-user. Developers of new algorithms
   may use it, but there is no guarantee that it will exist in the future.
 }

 \usage{\method{fitMultiscanAffine}{RGData}(this, slides=NULL, channels=c("R", "G"), aShift=NULL, method="L1", maxIter=30, reps=0.02, acc=1e-04, satSignal=2^16 - 1, bootstrap=FALSE, R=NULL, verbose=FALSE, ...)}

 \arguments{
  \item{slides}{Slides to be used in the fit \emph{and} that are to be
   calibrated. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{channels}{Specifies which channels to be fitted. If \code{\link[base]{NULL}}, all
   channels are considered.}
  \item{method}{\code{\link[base]{character}} specifying method for IWPCA fitting.
   If \code{"L1"}, the distances are minimised in \eqn{L_1}.
   No other methods are available.}
  \item{maxIter}{The maximum number of re-weigthed iterations.}
  \item{acc}{The accuracy for the stopping criterion.}
  \item{satSignal}{Signals equal to or above this threshold will not
    be used in the fitting.}
  \item{bootstrap}{If \code{\link[base:logical]{TRUE}}, bootstrap is applied to estimate the
    mean and standard deviation of the bootstraped estimates.}
  \item{R}{The number of bootstrap samples. Only effective if
   \code{bootstrap==TRUE}.}
  \item{reps}{Small value, which is added to the absolute value of
    the residuals to avoid "1/0" weights, i.e. "(1/0+reps)".}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, extra information is printed while running.}
  \item{...}{Additional arguments accepted by \code{\link[aroma.light]{iwpca}}.
   For instance, a N \code{\link[base]{vector}} of weights for each observation may be
    given, otherwise they get the same weight.}
 }

 \value{
   Returns a \code{\link[base]{list}} of fitted \code{\link{AffineModelFit}} objects for
   each channel.
 }

 \details{
  Fitting is done by iterated re-weighted principal component analysis
  (IWPCA).
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
   [1] H. Bengtsson, J. Vallon-Christersson and G. Jönsson, \emph{Calibration and assessment of channel-specific biases in microarray data with extended dynamical range}, BMC Bioinformatics, 5:177, 2004.
 \cr
 }

 \seealso{
   To calibrate the data using these fits see \code{\link[aroma:calibrateMultiscan.RGData]{*calibrateMultiscan}()}.
   For a similar strategy for between channel normalization see
   \code{\link[aroma:normalizeAffine.RGData]{*normalizeAffine}()}.
   For more information see \code{\link{RGData}}.
 }
\keyword{internal}
\keyword{methods}
