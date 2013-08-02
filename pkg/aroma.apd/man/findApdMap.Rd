%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findApdMap.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findApdMap}
\alias{findApdMap.default}
\alias{findApdMap}

\title{Search for APD map files in multiple directories}

\description{
 Search for APD map files in multiple directories.
}

\usage{
\method{findApdMap}{default}(mapType=NULL, paths=NULL, pattern="[.](a|A)(p|P)(m|M)$", ...)
}

\arguments{
 \item{mapType}{A \code{\link[base]{character}} string of the map type to search for.}
 \item{paths}{A \code{\link[base]{character}} \code{\link[base]{vector}} of paths to be searched.
   The current directory is always searched at the beginning.
   If \code{\link[base]{NULL}}, default paths are searched.  For more details, see below.}
 \item{pattern}{A regular expression file name pattern to match.}
 \item{...}{Additional arguments passed to \code{\link[R.utils]{findFiles}}.}
}

\value{
 Returns a \code{\link[base]{vector}} of the full pathnames of the files found.
}

\details{
  Note, the current directory is always searched at the beginning.
  This provides an easy way to override other files in the search path.
  If \code{paths} is \code{\link[base]{NULL}}, then a set of default paths are searched.
  The default search path is consituted of:
  \enumerate{
   \item \code{"."}
   \item \code{getOption("AFFX_APD_PATH")}
   \item \code{Sys.getenv("AFFX_APD_PATH")}
  }

  One of the easiest ways to set system variables for \R is to
  set them in an \code{.Renviron} file, see \code{\link[base]{Startup}}
  for more details.
}

\author{Henrik Bengtsson}



\keyword{file}
\keyword{IO}
