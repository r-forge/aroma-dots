%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SpotData.R
% 
% on Tue Jan 15 18:37:20 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getForegroundSD.SpotData}
\alias{getForegroundSD.SpotData}
\alias{SpotData.getForegroundSD}
\alias{getForegroundSD.SpotData}
\alias{getForegroundSD,SpotData-method}


\title{Gets (an approximation of) the standard deviation of the foreground pixels}

\usage{\method{getForegroundSD}{SpotData}(this, ...)}

\description{
 Gets (an approximation of) the standard deviation of the foreground pixels.
 Since the Spot software does not return estimates of the standard
 deviation of the foreground pixels, but the interquartile range (IQR) of
 their logarithm instead, the value returned by this method is only an
 approximation of the standard deviation.
 \emph{It is based on the assumption that the noise is symmetric on the
 non-logaritmic scale}.
}

\value{
  Returns a list of matrices that contain approximations of the estimated
  standard deviation of the pixels in the foreground region of the spots.
}

\details{
 By assuming that the noise is symmetric on the non-logaritmic scale we
 can first derive the IQR on the non-logaritmic scale since we know the
 median pixel value too. Let the non-logaritmic pixel value be \eqn{X} and
 the let \eqn{x=log[2](X)}. More over, let \eqn{X[0.50]} be the median
 value of all \eqn{X}'s, \eqn{X[0.25]} and \eqn{X[0.75]} be the 0.25 and
 0.75 quantile, respectively. \eqn{x[0.25], x[0.50], x[0.75]} are defined
 analogously. Assuming symmetric noise on the non-log scale means that
 \eqn{X[0.50]-X[0.25] == X[0.75]-X[0.50]} is true.
 Let \eqn{dX = X[0.75]-X[0.50]}. Given the IQR or the log pixel values,
 i.e. \eqn{xIQR=x[0.75]-x[0.25]}, one can derive
 \deqn{
   dX = \frac{2^{x_{IQR}}-1}{2^{x_{IQR}}+1} \cdot X_{0.50}
 }{
   dX = (2^xIQR-1/(2^xIQR+1)*X[0.50]
 }
 Hence, the IQR of the non-logarithmic pixel values can then be
 approximated by \eqn{2*dX}. If the noise is not only symmetric,
 but also Gaussian, the standard deviation is \eqn{1.35*2*dX},
 which is the value returned.
}

\examples{
  spot <- SpotData$read("spot123.spot", path=system.file("data-ex", package="aroma"))
  raw <- getRawData(spot)
  sd <- getForegroundSD(spot)
  raw$RSD <- sd$RSD; raw$GSD <- sd$GSD;
  subplots(4)
  plot(raw, "RSDvsR", col="red")
  plot(raw, "GSDvsG", col="green")
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{SpotData}}.
}
\keyword{internal}
\keyword{methods}
