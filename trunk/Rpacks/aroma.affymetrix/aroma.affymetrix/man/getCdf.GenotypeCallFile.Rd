%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getCdf.GenotypeCallFile}
\alias{getCdf.GenotypeCallFile}
\alias{GenotypeCallFile.getCdf}
\alias{getCdf.GenotypeCallFile}
\alias{getCdf,GenotypeCallFile-method}


\title{Gets the CDF}

\description{
 Gets the CDF for this call file.
}

\usage{\method{getCdf}{GenotypeCallFile}(this, ...)}

\arguments{
  \item{...}{Not used.}
}

\value{
 Returns an \code{\link{AffymetrixCdfFile}}.
}

\details{
 The CDF is looked up from the chip type attribute part of the file
 header.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{GenotypeCallFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
