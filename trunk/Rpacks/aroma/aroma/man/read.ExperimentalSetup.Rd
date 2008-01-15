%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ExperimentalSetup.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.ExperimentalSetup}
\alias{read.ExperimentalSetup}
\alias{ExperimentalSetup.read}
\alias{read.ExperimentalSetup}
\alias{read,ExperimentalSetup-method}


\title{Reads a file specifying the setup of a microarray experiment}

\description{
 Static method that read a experiment setup file describing a cDNA
 microarray experiment.
}

\synopsis{read.ExperimentalSetup(static, filename, path=NULL, ...)}
\usage{ExperimentalSetup$read(filename, path=NULL, ...)}

\arguments{
  \item{filename}{The filename of the setup file to be read.}
  \item{path}{Optional path to the setup file.}
  \item{verbose}{If \\code{\link[base:logical]{TRUE}}, information will printed out during
                 the file reading/parsing.}
}

\value{
  Returns an \code{\link{ExperimentalSetup}} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   setup <- ExperimentalSetup$read("ScanAlyze.setup", path=system.file("misc", package="aroma"))
   print(setup)
}

\seealso{
  For more information see \code{\link{ExperimentalSetup}}.
}
\keyword{internal}
\keyword{methods}
