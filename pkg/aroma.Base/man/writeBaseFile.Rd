%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  writeBaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{writeBaseFile}
\alias{writeBaseFile.default}
\alias{writeBaseFile}

\title{Low-level function to write a BASE file structure to a connection or a file}

\description{
 Low-level function to write a BASE file structure to a connection or a file.  Note that it is required that the structure is correct;
 Only minimal validation of structure is done, e.g. the existance of
 a mandatory header named 'section' is asserted.
}

\usage{\method{writeBaseFile}{default}(con, base, verbose=FALSE, ...)}

\arguments{
  \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} string filename.}
  \item{base}{The BASE file \code{\link[base]{list}} structure to be written.}
  \item{verbose}{Either a \code{\link[base]{logical}}, a \code{\link[base]{numeric}}, or a \code{\link[R.utils]{Verbose}}
    object specifying how much verbose/debug information is written to
    standard output. If a Verbose object, how detailed the information is
    is specified by the threshold level of the object. If a numeric, the
    value is used to set the threshold of a new Verbose object. If \code{\link[base:logical]{TRUE}},
    the threshold is set to -1 (minimal). If \code{\link[base:logical]{FALSE}}, no output is written
    (and neither is the \link[R.utils:R.utils-package]{R.utils} package required).
  }
  \item{...}{Not used.}
}

\value{
  Returns \code{\link[base:logical]{TRUE}}, if structure was succesfully written, otherwise \code{\link[base:logical]{FALSE}}.
}

\seealso{
  \code{\link{readBaseFile}}().
  \code{\link{writeBaseFileSection}}().
  See \code{stdout()} in \link[base]{showConnections} to write to standard
  output.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}



\keyword{file}
\keyword{IO}
