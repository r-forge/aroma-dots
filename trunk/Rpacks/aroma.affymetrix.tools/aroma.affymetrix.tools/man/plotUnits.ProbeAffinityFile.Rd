\name{plotUnits.ProbeAffinityFile}
\alias{plotUnits.ProbeAffinityFile}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting Function for ProbeAffinityFile }
\description{
    Plots probe affinities, usually only called to plot on top of existing function or calculate certain terms of model
}
\usage{
plotUnits.ProbeAffinityFile(this, units, intercept = 0, x.pos = NULL, add = T, type = "l", joinGroups = T, plot.it = T, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ ProbeAffinityFile object}
  \item{units}{ unit to plot }
  \item{intercept}{ Give intercept term (see details)}
  \item{x.pos}{ positions at which to plot the probe affinities, defaults to integers 1 to number of probes}
  \item{add}{ whether to add on top of existing plot }
  \item{type}{ type of plotting }
  \item{joinGroups}{ whether to join groups together }
  \item{plot.it}{ whether to plot}
  \item{\dots}{ additional arguments passed to \code{plot} or \code{lines}, depending on \code{add}. }
}
\details{
    This function is intended to allow the plotting of probe affinities atop the original expression values. Therefore, the option \code{intercept} allows for giving an amount to add to all the probe affinities so that they are on the right scale for the original intensities (the mean chip effect, usually)
}
\value{
  returns vector of probe affinities
}
\author{ Elizabeth Purdom}
