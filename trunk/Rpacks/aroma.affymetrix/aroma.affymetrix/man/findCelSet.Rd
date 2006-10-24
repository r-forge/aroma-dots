%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findCelSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findCelSet}
\alias{findCelSet.default}
\alias{findCelSet}

\title{Searches for one or several datasets}

\description{
  Searches for one or several datasets that matches certain criterias.
}

\usage{\method{findCelSet}{default}(name=NULL, chipType=NULL, paths=NULL, minCount=1, ..., firstOnly=TRUE, verbose=FALSE)}

\arguments{
 \item{name}{A regular expression pattern to match the name of the
   \emph{data set} subdirectory.}
 \item{chipType}{A regular expression pattern to match the name of
   the \emph{chip type} subdirectory.}
 \item{paths}{A \code{\link[base]{character}} \code{\link[base]{vector}} of paths to be searched.
   If \code{\link[base]{NULL}} certain default directories are scanned. See below.}
 \item{minCount}{The minimum number of CEL files in a directory for
   it to be considered a data set.}
 \item{...}{Not used.}
 \item{firstOnly}{If \code{\link[base:logical]{TRUE}}, the first possible match is returned,
   otherwise all matches are returned.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns \code{\link[base]{character}} \code{\link[base]{vector}} of full pathname to matching dataset.
}

\details{
  This function searches for data sets with a pathname pattern
  \preformatted{
     path/to/<data set name>/chip_data/<chip type>/
  }
  of a directory containing CEL files.
}

\section{Default search paths}{
  Note, the current directory is always searched at the beginning.
  This provides an easy way to override other files in the search path.
  If \code{paths} is \code{\link[base]{NULL}}, then a set of default paths are searched.
  The default search path is consituted of:
  \enumerate{
   \item \code{"."}
   \item \code{getOption("AFFX_DATA_PATH")}
   \item \code{Sys.getenv("AFFX_DATA_PATH")}
   \item \code{"data/"}
  }

  One of the easiest ways to set system variables for \R is to
  set them in an \code{.Renviron} file, e.g.
  \preformatted{
    # affxparser: Set default CDF path
    AFFX_DATA_PATH=${AFFX_DATA_PATH};M:/data
  }
  See \code{\link[base]{Startup}} for more details.
}

\examples{
# Search for a specific dataset
path <- findCelSet("SlaterH_etal_2004")
print(path)

# Search for a specific dataset of a specific chip type
path <- findCelSet("SlaterH_etal_2004", chipType="Xba")
print(path)

# Find a dataset with a specific chip type
path <- findCelSet(chipType="Sty", firstOnly=TRUE)
print(path)

# Find a dataset with at least 8 arrays
path <- findCelSet(minCount=8, firstOnly=TRUE)
print(path)

}

\author{Henrik Bengtsson \url{http://www.braju.com/R/}}



\keyword{file}
\keyword{IO}
