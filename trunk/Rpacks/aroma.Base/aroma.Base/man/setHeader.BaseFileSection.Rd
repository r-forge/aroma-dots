%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{setHeader.BaseFileSection}
\alias{setHeader.BaseFileSection}
\alias{BaseFileSection.setHeader}
\alias{setHeader.BaseFileSection}
\alias{setHeader,BaseFileSection-method}


 \title{Sets a header of a BASE section}

 \description{
   Sets a header of a BASE section.
 }

 \usage{\method{setHeader}{BaseFileSection}(this, name, value, ...)}

 \arguments{
  \item{name}{A \code{\link[base]{character}} string of the name of the header.}
  \item{value}{The value of the header. May be a \code{\link[base]{vector}}.}
  \item{type}{A \code{\link[base]{character}} string specifying the new type.}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the old header value, or \code{\link[base]{NULL}} if missing.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
