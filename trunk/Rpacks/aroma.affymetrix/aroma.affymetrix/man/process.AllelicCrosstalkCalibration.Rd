%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AllelicCrosstalkCalibration.R
% 
% on Sat Jan 06 01:48:24 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{process.AllelicCrosstalkCalibration}
\alias{process.AllelicCrosstalkCalibration}
\alias{AllelicCrosstalkCalibration.process}
\alias{process.AllelicCrosstalkCalibration}
\alias{process,AllelicCrosstalkCalibration-method}


\title{Calibrates the data set}

\description{
 Calibrates the data set.
}

\usage{\method{process}{AllelicCrosstalkCalibration}(this, ..., force=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Not used.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, data already calibrated is re-calibrated,
      otherwise not.}
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
# Define an example dataset
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Find any dataset
path <- findCelSet(minCount=5)
if (is.null(path))
  break

ds <- AffymetrixCelSet$fromFiles(path)
print(ds)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsQN <- normalizeQuantile(ds, subsetToAvg=0.50,
                               typesToUpdate="pm", verbose=verbose)
print(dsQN)

} # for (zzz in 0)
rm(zzz)
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AllelicCrosstalkCalibration}}.
}
\keyword{internal}
\keyword{methods}
