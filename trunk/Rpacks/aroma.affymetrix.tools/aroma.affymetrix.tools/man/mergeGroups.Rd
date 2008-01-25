\name{mergeGroups}
\alias{mergeGroups.AffymetrixCelFile}
\alias{mergeGroups.list}
\alias{mergeGroups.ProbeLevelModel}
\alias{mergeGroups.AffymetrixCelSet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
\method{mergeGroups}{AffymetrixCelFile}(this, units, fields = NULL, force = FALSE, logscale = FALSE, ...)
\method{mergeGroups}{AffymetrixCelSet}(this, units, arrays = NULL, fields = NULL, force = FALSE, logscale = FALSE, ...)
\method{mergeGroups}{ProbeLevelModel}(this, units, arrays = NULL, fields = NULL, force = FALSE, logscale = FALSE, ...)
\method{mergeGroups}{list}(this, arrays = NULL, fields = NULL, logscale = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this}{ object}
  \item{units}{ units which you want to merge groups }
  \item{fields}{ field to get to merge (can be more than one); default is 'intensities' }
  \item{arrays}{ indices of which arrays to included}
  \item{force}{ passed to \code{readUnits} commands}
  \item{logscale}{ logical; whether to take log2 of data }
  \item{\dots}{parameters passed to \code{getUnitMergeGroupsFunction}; not used currently}
}
\details{
    If \code{this} is a list, must be in same structure as returned from \code{readUnits}.
}
\value{
    Returns list in same structure as \code{readUnits} only with groups merged (by rows)
}
\author{ Elizabeth Purdom}
