\name{plotUnits.RmaPlm}
\alias{plotUnits.RmaPlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot the data and estimated effects from RmaPlm fit}
\description{
    Plot (using matplot) the original (log2) intensities and various aspects of the fitted model.
}
\usage{
plotUnits.RmaPlm(this, units, arrays = NULL, joinGroups = T, includeProbeEffect = F, includeChipEffect = T, plotRaw = F, plotProbeEffect = F, probePlotArgs = list(lwd = 3, lty = 1, col = "black"), plotByArray = T, plotIntercept = F, interceptPlotArgs = list(col = "grey", lwd = 2, lty = 2), plot.it = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ RmaPlm object }
  \item{units}{ unit to plot (must be just 1) }
  \item{arrays}{ vector of indices of which arrays to plot (in order of the arrays in this)}
  \item{joinGroups}{ logical; whether to connect the lines that cross group boundaries}
  \item{includeProbeEffect}{ logical; if T, plots original intensities. If F, plots the intensities with probe effect subtracted off}
  \item{includeChipEffect}{ logical; if T, plots original intensities. If F, plots the intensities with chip effect subtracted off}
  \item{plotRaw}{ logical; whether to plot raw data (over rides other options) }
  \item{plotProbeEffect}{ logical; whether to superimpose the probe effect (affinity) estimates }
  \item{probePlotArgs}{ list of named arguments to pass to plotUnits.ProbeAffinityFile (if plotProbeEffect=T) }
  \item{plotByArray}{ logical; if true, each probe is separate line and axis is samples (default is other way) }
  \item{plotIntercept}{ logical; whether to superimpose the intercept (mean of chip effect) }
  \item{interceptPlotArgs}{ list of named arguments to pass to \code{abline} for plotting intercept}
  \item{plot.it}{ whether to plot }
  \item{\dots}{ arguments to pass to plotUnits.AffymetrixCelSet}
}
\value{
  returns matrix of data plotted (with effects subtracted off as requested)
}
\author{ Elizabeth Purdom}
