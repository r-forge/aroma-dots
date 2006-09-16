%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{plotMvsPosition.AffymetrixCelFile}
\alias{plotMvsPosition.AffymetrixCelFile}
\alias{AffymetrixCelFile.plotMvsPosition}
\alias{plotMvsPosition.AffymetrixCelFile}
\alias{plotMvsPosition,AffymetrixCelFile-method}


\title{Plots log-ratio versus physical position in a scatter plot}

\description{
 Plots log-ratio versus physical position in a scatter plot.
}

\usage{\method{plotMvsPosition}{AffymetrixCelFile}(this, reference, chromosome, region=NULL, gdas, ylab=ylab, xlab="Physical position [Mb]", xlim=NULL, ..., what="M", annotate=TRUE)}

\arguments{
  \item{reference}{The reference channel, i.e. the denominator in the
    log ratios.}
  \item{x}{The other variable.  A \code{\link[base]{double}} \code{\link[base]{vector}}.}
  \item{chromosome}{The chromosome to be plotted.}
  \item{region}{The region to be plotted.  If \code{\link[base]{NULL}}, the whole chromosome
    is plotted.}
  \item{pch}{The plot symbol.}
  \item{xlim,ylim}{The range of the x and the y axes.}
  \item{xlab,ylab}{The labels on the x and the y axes.}
  \item{...}{Additional arguments passed to \code{\link[graphics]{plot}}.}
  \item{annotate}{If \code{\link[base:logical]{TRUE}}, the plot is annotated with information about
    the data plotted, otherwise not.}
}

\value{
 Returns (invisibly) a \code{\link[base]{data.frame}} with the plotted data in the
 first two columns, and remaining data in the following columns.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:plotMvsX.AffymetrixCelFile]{*plotMvsX}()}.
  \code{\link[aroma.affymetrix:plotMvsA.AffymetrixCelFile]{*plotMvsA}()}.
  \code{\link[aroma.affymetrix:smoothScatterMvsA.AffymetrixCelFile]{*smoothScatterMvsA}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}
\keyword{internal}
\keyword{methods}
