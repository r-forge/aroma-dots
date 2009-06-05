%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getSections.BaseFile}
\alias{getSections.BaseFile}
\alias{BaseFile.getSections}
\alias{getSections.BaseFile}
\alias{getSections,BaseFile-method}


 \title{Gets a subset or all BASE sections}

 \description{
   Gets a subset or all BASE sections specified by their names and/or their header contents.
 }

 \usage{\method{getSections}{BaseFile}(this, types=NULL, headerValues=NULL, regexpr=FALSE, ...)}

 \arguments{
  \item{types}{A \code{\link[base]{vector}} of \code{\link[base]{character}} strings specifying which sections
    to retrieve.}
  \item{headers}{A \code{\link[base]{vector}} of \code{\link[base]{character}} strings specifying headers that
    the sections must contain in order to be returned.}
  \item{regexpr}{If \code{\link[base:logical]{TRUE}}, the \code{types} and \code{headers} are
    interpreted as regular expressions, other exact matching is required.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{list}} of \code{\link{BaseFileSection}}s.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
