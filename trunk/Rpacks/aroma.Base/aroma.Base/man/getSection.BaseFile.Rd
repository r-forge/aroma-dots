%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{getSection.BaseFile}
\alias{getSection.BaseFile}
\alias{BaseFile.getSection}
\alias{getSection.BaseFile}
\alias{getSection,BaseFile-method}


 \title{Gets one section of this BASE file structure}

 \description{
   Gets one section of this BASE file structure.  If more than one exists of the same type, only the first
   is returned.
 }

 \usage{\method{getSection}{BaseFile}(this, type, ...)}

 \arguments{
  \item{type}{A \code{\link[base]{character}} string specifying the section to retrieve.}
  \item{...}{Arguments passed to \code{\link[aroma.Base:getSections.BaseFile]{*getSections}()}.}
 }

 \value{
   Returns a \code{\link{BaseFileSection}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
