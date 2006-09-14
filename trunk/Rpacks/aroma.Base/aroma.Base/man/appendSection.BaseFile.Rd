%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{appendSection.BaseFile}
\alias{appendSection.BaseFile}
\alias{BaseFile.appendSection}
\alias{appendSection.BaseFile}
\alias{appendSection,BaseFile-method}


 \title{Adds a section to this BASE file structure}

 \description{
   Adds a section to this BASE file structure.
 }

 \usage{\method{appendSection}{BaseFile}(this, section, ...)}

 \arguments{
  \item{section}{A \code{\link{BaseFileSection}} to be appended.}
  \item{...}{Not used.}
 }

 \value{
   Returns \code{\link[base:logical]{TRUE}} if section was added.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma.Base:getSection.BaseFile]{*getSection}()}.
   \code{\link[aroma.Base:removeSection.BaseFile]{*removeSection}()}.
   \code{\link[aroma.Base:replaceSection.BaseFile]{*replaceSection}()}.
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
