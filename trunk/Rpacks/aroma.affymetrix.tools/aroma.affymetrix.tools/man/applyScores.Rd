\name{applyScores}
\alias{applyScores}
\alias{flatScores}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Apply for Saved Scores }
\description{
Apply a function to the group values returned from calculating scores}
\usage{
applyScores(filename, path = NULL, FUN = function(x) {x})
}
\usage{
flatScores(filename, path = NULL, FUN = function(x) {
    x
}, componentNames = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{ the core part of filename of saved files (given to calculateScores) }
  \item{path}{ (optional) the path for the folder where the scores are saved }
  \item{FUN}{ a function that acts on the group elements of the scores list }
  \item{componentNames}{character string giving the element of the named group list to be extracted before flattening scores}

}
\details{
applyScores applies FUN to each group element (so FUN needs to take into account sublist structure, if any, in the scores objects).
Returns with unit/group structure, but perhaps not further sublist (as would be in AffymetrixCelSet) depending on FUN. Default just recreates entire list object.

flatScores returns all of the scores without list structure. You can give the option as to which element of the named group list should be kept, in case there are different categories of the scores saved in the list. FUN is then applied on the entire flatten list, so must work on vectors. Similar output could be achieved with the appropriate FUN in applyScores.
}
\value{
    List in AffymetrixCelSet structure (but see details)
}
\author{Elizabeth Purdom}
