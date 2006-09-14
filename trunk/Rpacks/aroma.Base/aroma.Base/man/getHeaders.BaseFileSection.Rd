%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getHeaders.BaseFileSection}
\alias{getHeaders.BaseFileSection}
\alias{BaseFileSection.getHeaders}
\alias{getHeaders.BaseFileSection}
\alias{getHeaders,BaseFileSection-method}


 \title{Gets the headers of a BASE section}

 \description{
   Gets the headers of a BASE section.
 }

 \usage{\method{getHeaders}{BaseFileSection}(this, names=NULL, regexpr=FALSE, ...)}

 \arguments{
  \item{names}{A \code{\link[base]{vector}} of \code{\link[base]{character}} strings specifying which headers
    to retrieve.}
  \item{regexpr}{If \code{\link[base:logical]{TRUE}}, the \code{names} are interpreted as regular
    expressions, other exact matching is required.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{list}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getHeader.BaseFileSection]{*getHeader}()}.
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
