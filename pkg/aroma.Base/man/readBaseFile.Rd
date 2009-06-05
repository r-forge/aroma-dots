%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readBaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readBaseFile}
\alias{readBaseFile.default}
\alias{readBaseFile}

\title{Low-level function to read a BASE file (stream) from a connection or a file}

\description{
 Low-level function to read a BASE file (stream) from a connection or a file.
}

\usage{\method{readBaseFile}{default}(con, ..., verbose=FALSE)}

\arguments{
  \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} string filename.}
  \item{...}{Arguments such as the important \code{extractSpotsData},
    passed to \code{\link{readBaseFileSection}}().}
  \item{verbose}{Either a \code{\link[base]{logical}}, a \code{\link[base]{numeric}}, or a \code{\link[R.utils]{Verbose}}
    object specifying how much verbose/debug information is written to
    standard output. If a Verbose object, how detailed the information is
    is specified by the threshold level of the object. If a numeric, the
    value is used to set the threshold of a new Verbose object. If \code{\link[base:logical]{TRUE}},
    the threshold is set to -1 (minimal). If \code{\link[base:logical]{FALSE}}, no output is written
    (and neither is the \link[R.utils:R.utils-package]{R.utils} package required).
  }
}

\value{
  Returns a named \code{\link[base]{list}} structure containing the parsed BASE structure.
  The names of the elements in the list are the same as the names (types)
  of the sections as given by section header 'section'.
}

\examples{
# Read a simple BASE file structure (adopted from [2])
filename <- system.file("data-ex", "baseExample.base", package="aroma.Base")
base <- readBaseFile(filename)
print(base)

# Write it back to file
tmpfile <- tempfile()
writeBaseFile(tmpfile, base)

# Read the newly create file
base2 <- readBaseFile(tmpfile)
file.remove(tmpfile)

# Compare with the original structure
stopifnot(equals(base, base2))
}

\seealso{
  \code{\link{writeBaseFile}}().
  \code{\link{readBaseFileSection}}().
  See \code{stdin()} in \link[base]{showConnections} to read from standard input.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  [1] BASE - BioArray Software Environment,\cr
      \url{http://base.thep.lu.se/}\cr

  [2] Carl Troein, How to write a plugin for BASE, April 2003.\cr
      \url{http://base.thep.lu.se/documentation/development/plugins.txt}\cr
      \url{http://opensource.microarray.omrf.org/wiki/bin/view/BASE/PluginWritingHowto}
}



\keyword{file}
\keyword{IO}
