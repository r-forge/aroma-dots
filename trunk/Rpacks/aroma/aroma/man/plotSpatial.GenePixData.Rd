%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.GenePixData}
\alias{plotSpatial.GenePixData}
\alias{GenePixData.plotSpatial}
\alias{plotSpatial.GenePixData}
\alias{plotSpatial,GenePixData-method}


\title{Creates a spatial plot of a slide}

\description{
 Creates a spatial plot of a slide.
}

\usage{\method{plotSpatial}{GenePixData}(this, what=NULL, slide=1, include=NULL, exclude=NULL, pch=20, yaxt=NULL, col="auto", palette="redgreen", A.range=c(0, 16), M.range=c(-1, 1), look=c("real", "classic"), cex=NULL, style=NULL, ...)}

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
  \item{yaxt, ...}{Parameters as accepted by \code{plot}.}
}

\value{Returns nothing.}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  subplots(2)

  opar <- par(bg="black")
  plotSpatial(gpr)
  par(opar)

  opar <- par(bg="black")
  plotSpatial(gpr, palette="blueyellow")
  par(opar)
}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
