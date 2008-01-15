%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GalLayout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.GalLayout}
\alias{read.GalLayout}
\alias{GalLayout.read}
\alias{read.GalLayout}
\alias{read,GalLayout-method}


\title{Reads a GenePix Array List (GAL) file}

\description{
 Reads a GenePix Array List (GAL) file.
}

\synopsis{read.GalLayout(this, filename, path=NULL, solve=FALSE, verbose=FALSE, ...)}
\usage{GalLayout$read(filename, path=NULL, solve=FALSE, verbose=FALSE, ...)}

\arguments{
  \item{filename}{The name of the file.}
  \item{path}{Optional path where the data should be written.}
  \item{solve}{If \code{\link[base:logical]{TRUE}}, the reading is slower but more forgiving to
        errors in the input file.}
  \item{verbose}{If \code{\link[base:logical]{TRUE}}, helpful information is printed at each step
        of the parsing of the input file.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{\dontrun{See help(GalLayout) for an example.}}

\seealso{
  For writing a GalLayout object back to a file \code{\link[aroma:write.GalLayout]{*write}()}.
  For more information see \code{\link{GalLayout}}.
}
\keyword{internal}
\keyword{methods}
