%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  parseTcgaDccArchiveName.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{parseTcgaDccArchiveName.character}
\alias{parseTcgaDccArchiveName.character}
\alias{parseTcgaDccArchiveName.character}


\title{Parses a TCGA DCC archive name}

\description{
  Parses a TCGA DCC archive name into its components.
}

\usage{\method{parseTcgaDccArchiveName}{character}(name, unlist=FALSE, ...)}

\arguments{
 \item{name}{A \code{\link[base]{character}} string.}
 \item{unlist}{A \code{\link[base]{logical}}.}
 \item{...}{Not used.}
}

\value{
  Returns a named \code{\link[base]{character}} \code{\link[base]{vector}} of length N=6.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  name <- "broad.mit.edu_GBM.HG-U133_Plus_2.1.0.0"
  parts <- parseTcgaDccArchiveName(name)
}

\references{
  [1] The TCGA Network, \emph{TCGA Data Primer - Version 1.0}, July 2008.\cr
}


\keyword{methods}
\keyword{file}
