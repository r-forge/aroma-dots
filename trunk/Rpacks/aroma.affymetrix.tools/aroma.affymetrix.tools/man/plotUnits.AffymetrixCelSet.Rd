\name{plotUnits.AffymetrixCelSet}
\alias{plotUnits.AffymetrixCelSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to plot AffymetrixCelSet data}
\description{
  Plots (using matplot) the data from a single unit
}
\usage{
\method{plotUnits}{AffymetrixCelSet}(this, units, dataMat = NULL, arrays = NULL, field = "intensities", markGroups = TRUE, markcol = "grey", labelGroups = FALSE, labelcol = "black", plotByArray = TRUE, marklwd = 1.5, joinGroups = TRUE, plot.it = TRUE, main = NULL, type = "l", xaxt = "n", ylab = NULL, xlab = "", logscale = TRUE, col = 1:6, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ AffymetixCelSet object}
  \item{units}{ unit to plot (must be just 1) }
  \item{dataMat}{ (optional) data matrix to plot instead of default matrix}
  \item{arrays}{ vector of indices of which arrays to plot (in order of the arrays in this)}
  \item{field}{ which field to plot (default 'intensities')}
  \item{markGroups}{ logical; whether to mark the groups with vertical lines (for exons)}
  \item{markcol}{ single color for the vertical lines marking groups}
  \item{marklwd}{ the width of the  vertical group markings}
  \item{labelGroups}{ label the horizontal axis with the group names}
  \item{labelcol}{ the colors for the group labels (can be vector)}
  \item{plotByArray}{ logical; if true, each probe is separate line and axis is samples (default is other way) }
  \item{joinGroups}{ whether to continue the lines between groups or keep them separate }
  \item{plot.it}{ whether to plot}
  \item{logscale}{logical; whether to take log2 of the data}
  \item{col}{ colors for each line (can be vector)}
  \item{xaxt,xlab,ylab,main,type}{ standard plotting options, though with some nonstandard defaults}
  \item{\dots}{ additional arguments passed to matplot}
}

\value{
  returns (invisibly) the data matrix from the AffymetrixCelSet that was plotted
}
\author{ Elizabeth Purdom}
