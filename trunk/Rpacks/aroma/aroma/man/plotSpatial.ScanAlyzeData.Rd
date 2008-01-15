%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyzeData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.ScanAlyzeData}
\alias{plotSpatial.ScanAlyzeData}
\alias{ScanAlyzeData.plotSpatial}
\alias{plotSpatial.ScanAlyzeData}
\alias{plotSpatial,ScanAlyzeData-method}


\title{Creates a spatial plot of a slide}

\description{
 Creates a spatial plot of a slide.
}

\usage{\method{plotSpatial}{ScanAlyzeData}(this, slide=1, pch=20, yaxt=NULL, col=NULL, palette="redgreen", A.range=c(0, 16), M.range=c(-1, 1), ...)}

\arguments{
  \item{slide}{The slide to be plotted.}
  \item{pch}{The spot symbol. Default is \code{20} (solid disk).}
  \item{col}{The color of the spots. If \code{\link[base]{NULL}}, default color is used.}
  \item{palette}{If \code{col} is not specified, colors are generated
  automaticially from the signals in the two foreground channels.
  If \code{redgreen}, a red to green colors scale will be used.
  If \code{blueyellow}, a blue to yellow colors scale will be used.}
  \item{A.range}{The range of the log (base 2) spot intensities. Used
   only for generating colors.}
  \item{M.range}{The range of the log (base 2) spot ratios. Used
   only for generating colors.}
  \item{yaxt, ...}{Parameters accepted by \code{plot}.}
}

\value{Returns nothing.}

\details{
  The ScanAlyze software does not return the center position of a spot,
  but the \code{TOP} \code{LEFT} \code{BOT} and \code{RIGHT} coordinates.
  This method assumes that the center of the spot is in the center of
  this box.
}

\examples{
  sa <- ScanAlyzeData$read("group4.dat", path=system.file("data-ex", package="aroma"))

  subplots(2)

  opar <- par(bg="black")
  plotSpatial(sa)
  par(opar)

  opar <- par(bg="black")
  plotSpatial(sa, palette="blueyellow")
  par(opar)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{ScanAlyzeData}}.
}
\keyword{internal}
\keyword{methods}
