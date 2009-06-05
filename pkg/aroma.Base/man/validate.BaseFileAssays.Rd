%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileAssays.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{validate.BaseFileAssays}
\alias{validate.BaseFileAssays}
\alias{BaseFileAssays.validate}
\alias{validate.BaseFileAssays}
\alias{validate,BaseFileAssays-method}


 \title{Validates assays section}

 \description{
   Validates assays section.
 }

 \usage{\method{validate}{BaseFileAssays}(this, action=c("error", "warning", "quiet"), ...)}

 \arguments{
  \item{action}{A \code{\link[base]{character}} string specifying action if a problem is
               detected. See \code{\link{validate.BaseFileSection}}() for details.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{character}} string explaining problems detected. If no problems
   where detected, \code{\link[base]{NULL}} is returned.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileAssays}}.
 }
\keyword{internal}
\keyword{methods}
