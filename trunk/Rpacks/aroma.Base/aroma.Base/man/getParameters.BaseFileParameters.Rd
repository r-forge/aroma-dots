%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileParameters.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getParameters.BaseFileParameters}
\alias{getParameters.BaseFileParameters}
\alias{BaseFileParameters.getParameters}
\alias{getParameters.BaseFileParameters}
\alias{getParameters,BaseFileParameters-method}


 \title{Gets parameters from a BASE parameters section}

 \description{
   Gets parameters from a BASE parameters section.  All parameters are retrieved via
   \code{\link[aroma.Base:getParameter.BaseFileParameters]{*getParameter}()} so that certain parameters are pre-processed.
 }

 \usage{\method{getParameters}{BaseFileParameters}(this, ...)}

 \arguments{
  \item{...}{Arguments passed to \code{\link[aroma.Base:getHeader.BaseFileParameters]{*getHeader}()}.}
 }

 \value{
   Returns a named \code{\link[base]{list}} structure of parameters.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getParameter.BaseFileParameters]{*getParameter}()}.
   For more information see \code{\link{BaseFileParameters}}.
 }
\keyword{internal}
\keyword{methods}
