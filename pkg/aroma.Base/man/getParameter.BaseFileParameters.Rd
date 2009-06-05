%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFileParameters.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getParameter.BaseFileParameters}
\alias{getParameter.BaseFileParameters}
\alias{BaseFileParameters.getParameter}
\alias{getParameter.BaseFileParameters}
\alias{getParameter,BaseFileParameters-method}


 \title{Gets one parameter from a BASE parameters section}

 \description{
   Gets one parameter from a BASE parameters section.
 }

 \usage{\method{getParameter}{BaseFileParameters}(this, name, ...)}

 \arguments{
  \item{name}{A \code{\link[base]{character}} string of the name of the parameter.}
  \item{...}{Not used.}
 }

 \value{
   Returns the value of a parameter.
 }

 \details{
   Certain parameters are pre-processed. This is a trial behavior. See
   source code, i.e. \code{print(getParameter.BaseFileParameters)}
   for details.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFileParameters}}.
 }
\keyword{internal}
\keyword{methods}
