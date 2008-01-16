%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.PLOT.R
% 
% on Tue Jan 15 18:33:31 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial3d.MicroarrayData}
\alias{plotSpatial3d.MicroarrayData}
\alias{MicroarrayData.plotSpatial3d}
\alias{plotSpatial3d.MicroarrayData}
\alias{plotSpatial3d,MicroarrayData-method}


\title{Plots a 3-dimensional spatial representation of a field}

\description{
 Plots a 3-dimensional spatial representation of a field.
}

\usage{\method{plotSpatial3d}{MicroarrayData}(this, field, slide=1, pch=176, phi=40, theta=-25, xlab="x", ylab="y", zlab=field, log=NULL, ...)}

\arguments{
 \item{field}{Name of the field to be plotted.}
 \item{slide}{Slide number to be plotted.}
 \item{pch}{Default value is 176 (small circle).}
 \item{theta, phi}{Angles defining the viewing direction.
       \code{theta} gives the azimuthal direction and
       \code{phi} the colatitude.}
 \item{xlab,ylab,zlab}{Labels for the x, the y, and the z axis.}
 \item{log}{The base of logarithm to be used for the z dimension.
            If \code{\link[base]{NULL}}, the logarithm is not calculated.}
 \item{...}{Other arguments accepted by \code{\link[R.basic]{plot3d}}.}
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)
  subplots(4)
  plotSpatial3d(raw, "Rb", pch=".", col="red")
  plotSpatial3d(raw, "Gb", pch=".", col="green")
  plotSpatial3d(ma, "M", pch=".")
  plotSpatial3d(ma, "A", pch=".")
}

\seealso{
  \code{\link[aroma:plotXY.MicroarrayData]{*plotXY}()} and \code{\link[aroma:plotSpatial.MicroarrayData]{*plotSpatial}()},
  which are called by \code{plot}.
  \code{\link[R.basic]{plot3d}}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
