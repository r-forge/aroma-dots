%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelSet.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{normalizeQuantile.AffymetrixCelSet}
\alias{normalizeQuantile.AffymetrixCelSet}
\alias{AffymetrixCelSet.normalizeQuantile}
\alias{normalizeQuantile.AffymetrixCelSet}
\alias{normalizeQuantile,AffymetrixCelSet-method}


\title{Normalizes samples to have the same empirical distribution}

\description{
 Normalizes samples to have the same empirical distribution.
}

\usage{\method{normalizeQuantile}{AffymetrixCelSet}(this, path=NULL, name="normQuantiles", subsetToAvg=subsetToUpdate, typesToAvg=typesToUpdate, subsetToUpdate=NULL, typesToUpdate=NULL, xTarget=NULL, ..., verbose=FALSE)}

\arguments{
  \item{path}{The path where to save the normalized data files.
    If \code{\link[base]{NULL}}, a default dataset and chip-type specific path is used.}
  \item{name}{The name of the normalized data set, which will also be
    part of the default path.}
  \item{subsetToAvg}{The probes to calculate average empirical
    distribution over.  If a single \code{\link[base]{numeric}} in (0,1), then this
    fraction of all probes will be used.
    If \code{\link[base]{NULL}}, all probes are considered.}
  \item{typesToAvg}{Types of probes to be used when calculating the
    average empirical distribution.
    If \code{"pm"} and \code{"mm"} only perfect-match and mismatch
    probes are used, respectively. If \code{"pmmm"} both types are used.
  }
  \item{subsetToUpdate}{The probes to be updated.
    If \code{\link[base]{NULL}}, all probes are updated.}
  \item{typesToUpdate}{Types of probes to be updated.}
  \item{xTarget}{A \code{\link[base]{numeric}} \code{\link[base]{vector}}.  The empirical distribution
    to which all arrays should be normalized to.}
  \item{...}{Additional arguments passed to \code{normalizeQuantile()}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns a \code{\link[base]{double}} \code{\link[base]{vector}}.
}

\examples{
for (zzz in 0) {

# Setup verbose output
verbose <- Arguments$getVerbose(-2)
timestampOn(verbose)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("ds")) {
  # Search for a chip_data/ directory...
  dataPath <- filePath("chip_data", expandLinks="any")
  dirs <- list.files(path=dataPath, full.names=TRUE)
  ds <- NULL
  for (path in dirs) {
    if (isDirectory(path)) {
      # Try to define a dataset
      tryCatch({
        ds <- AffymetrixCelSet$fromFiles(path)
        verbose && cat(verbose, "Found a dataset: ", path)
        break
      }, error = function(ex) {})
    }
  }
  if (is.null(ds))
    break
}

verbose && cat(verbose, "Working with dataset:")
print(ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantiles(ds, subsetToAvg=0.50,
                               typesToUpdate="pm", verbose=verbose)
print(dsQN)


} # for (zzz in 0)
rm(zzz)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.light]{normalizeQuantile.numeric}}
  For more information see \code{\link{AffymetrixCelSet}}.
}
\keyword{internal}
\keyword{methods}
