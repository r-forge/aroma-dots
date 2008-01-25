\name{plotUnits.FirmaModel}
\alias{plotUnits.FirmaModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting function for Firma results }
\description{
    Plots image plot of the firma scores, per exon, as well as the original data and the weights, aligned by probe.
}
\usage{
\method{plotUnits}{FirmaModel}(this, units, arrays = NULL, colSchemeProbe = brewer.pal(11, "RdBu"), colSchemeFIRMA = c("white", brewer.pal(9, "YlGnBu")), breaksProbe = NULL, breaksFIRMA = NULL, colArrays = NULL, arrayFactor = NULL, arrayNames = NULL, mar = c(2.1, 7.1, 0.5, 0.2), intersp = 2.1, plotResiduals = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ FirmaModel object }
  \item{units}{ unit to plot (must be just 1) }
  \item{arrays}{ vector of indices of which arrays to plot (in order of the arrays in this)}
  \item{colSchemeProbe}{ color palette to use for the probe-level residuals image}
  \item{colSchemeFIRMA}{ color palette to use for the FIRMA scores }
  \item{breaksProbe}{ breaks for the probe-level residuals (passed to \code{image}) }
  \item{breaksFIRMA}{ breaks for the FIRMA scores (passed to \code{image}) }
  \item{colArrays}{ color of the arrays }
  \item{arrayFactor}{ factor for grouping the arrays in the image plots}
  \item{arrayNames}{ Names for the arrays}
  \item{mar}{ global \code{par} option for space on the edges of the set of aligned plots}
  \item{intersp}{ space to be put between the plots; passed to the \code{par} option \code{mar} }
  \item{plotResiduals}{ logical, whether to give an image plot of the probe level residuals }
  \item{\dots}{ arguments passed to \code{plotUnits.RmaPlm}}
}
\details{
Lines will be drawn in the image plots to separate the groups assigned by the factor variable \code{arrayFactor}. If \code{colArrays} is not given, \code{arrayFactor} will give color assignments in accordence with the factor grouping there. If \code{colArrays} is given, it will override this choice. This allows a separate grouping of the arrays different from the coloring.
}
\value{
Nothing is returned
}
\author{ Elizabeth Purdom}
