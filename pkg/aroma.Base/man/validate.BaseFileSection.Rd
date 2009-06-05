%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileSection.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{validate.BaseFileSection}
\alias{validate.BaseFileSection}
\alias{BaseFileSection.validate}
\alias{validate.BaseFileSection}
\alias{validate,BaseFileSection-method}


 \title{Validates BASE file section}

 \description{
   Validates BASE file section.
 }

 \usage{\method{validate}{BaseFileSection}(this, action=c("error", "warning", "none"), ...)}

 \arguments{
  \item{action}{A \code{\link[base]{character}} string specifying action if a problem is
               detected. If \code{"error"}, an exception is thrown,
               if \code{"warning"}, a warning is generated, and if
               \code{"none"}, to action is taken.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{character}} string explaining problems detected. If no problems
   where detected, \code{\link[base]{NULL}} is returned.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileSection}}.
 }
\keyword{internal}
\keyword{methods}
