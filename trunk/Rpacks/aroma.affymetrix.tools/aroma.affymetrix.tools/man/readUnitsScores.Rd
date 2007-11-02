\name{readUnitsScores}
\alias{readUnitsScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Read Scores for Units}
\description{
Like readUnits functions, finds scores for given unit names
}
\usage{
readUnitsScores(names, filename, path = NULL, verbose = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{names}{ names of the units to read the scores of (note different from readUnits functions)}
   \item{filename}{ the core part of filename of saved files (given to calculateScores) }
  \item{path}{ (optional) the path for the folder where the scores are saved }
}
\details{
}
\value{
    Returns the elements of the scores' list structure corresponding to unit names
}
\author{Elizabeth Purdom}
