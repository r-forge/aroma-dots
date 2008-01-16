%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QuantArrayData.R
% 
% on Tue Jan 15 18:35:00 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotSpatial.QuantArrayData}
\alias{plotSpatial.QuantArrayData}
\alias{QuantArrayData.plotSpatial}
\alias{plotSpatial.QuantArrayData}
\alias{plotSpatial,QuantArrayData-method}


\title{Creates a spatial plot of a slide}

\description{
 Creates a spatial plot of a slide.
}

\usage{\method{plotSpatial}{QuantArrayData}(this, what=NULL, slide=1, pch=20, yaxt=NULL, col="auto", palette="redgreen", A.range=c(0, 16), M.range=c(-1, 1), style=c("real", "classic"), ...)}

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
 \dontrun{
  # At the moment there is no QuantArray sample files in the package...
  # qa <- QuantArrayData$read("quantarray123.txt", path=system.file("data-ex", package="aroma"))
  # ...will use a GenePix sample file instead.
  qa <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  subplots(2)

  opar <- par(bg="black")
  plotSpatial(qa)
  par(opar)

  opar <- par(bg="black")
  plotSpatial(qa, palette="blueyellow")
  par(opar)
 }
}

\seealso{
  For more information see \code{\link{QuantArrayData}}.
}
\keyword{internal}
\keyword{methods}
