%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenotypeCallSet.R
% 
% on Sat Jan 06 01:52:11 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setCdf.GenotypeCallSet}
\alias{setCdf.GenotypeCallSet}
\alias{GenotypeCallSet.setCdf}
\alias{setCdf.GenotypeCallSet}
\alias{setCdf,GenotypeCallSet-method}


\title{Sets the CDF structure for this genotype call set}

\description{
 Sets the CDF structure for this genotype call set.  This structures is assigned to all genotype call files
   in the set.
}

\usage{\method{setCdf}{GenotypeCallSet}(this, cdf, ...)}

\arguments{
  \item{cdf}{An \code{\link{AffymetrixCdfFile}} object.}
  \item{...}{Not used.}
}

\value{
 Returns nothing.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:getCdf.GenotypeCallSet]{*getCdf}()}.
  For more information see \code{\link{GenotypeCallSet}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
