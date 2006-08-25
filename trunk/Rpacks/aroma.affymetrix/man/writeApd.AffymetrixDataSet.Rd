%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixDataSet.WRITE.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{writeApd.AffymetrixDataSet}
\alias{writeApd.AffymetrixDataSet}
\alias{AffymetrixDataSet.writeApd}
\alias{writeApd.AffymetrixDataSet}
\alias{writeApd,AffymetrixDataSet-method}


\title{Generates APD files for each of the data files in this data set}

\description{
  Generates APD files for each of the data files in this data set.
}

\usage{\method{writeApd}{AffymetrixDataSet}(this, apdMap="byChipType", ...)}

\arguments{
  \item{mapType}{The type of read map for the generated APD file.
    If \code{\link[base]{NULL}}, no remapping of the cell indices is done.
    If \code{"asChipType"}, the map type is the same as the chip type
    of the CEL file.
    If any other \code{\link[base]{character}} string, it sets the map type to that string.
    Note that there must be a APD map file with this type that can
    be found by \code{\link[aroma.apd]{findApdMap}}.
  }
 \item{...}{Not used.}
}

\value{
  Returns a new \code{AffymetrixDataSet} instance.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{AffymetrixDataSet}}.
}
\keyword{internal}
\keyword{methods}
