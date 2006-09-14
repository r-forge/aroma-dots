%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  extractBaseFileSpotsData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extractBaseFileSpotsData}
\alias{extractBaseFileSpotsData.default}
\alias{extractBaseFileSpotsData}

\title{Low-level function to extract assay data and save it to file instead}

\description{
 Low-level function to extract assay data and save it to file instead to save memory.
}

\usage{\method{extractBaseFileSpotsData}{default}(section, serialize=TRUE, ...)}

\arguments{
  \item{section}{A \code{\link[base]{list}} structure containing a 'spots' section.}
  \item{serialize}{If \code{\link[base:logical]{TRUE}}, non-serialized data is serialized before being
    extracted so that data is saved in seperate files for each assay.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} with element \code{section} containing the modified
  'spots' section.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link{readBaseFileSection}}().
}



\keyword{file}
\keyword{IO}
