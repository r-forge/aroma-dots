%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BaseFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{replaceSection.BaseFile}
\alias{replaceSection.BaseFile}
\alias{BaseFile.replaceSection}
\alias{replaceSection.BaseFile}
\alias{replaceSection,BaseFile-method}


 \title{Replaces a section in a BASE file structure}

 \description{
   Replaces a section in a BASE file structure.
 }

 \usage{\method{replaceSection}{BaseFile}(this, oldSection, newSection, append=FALSE, ...)}

 \arguments{
  \item{oldSection}{The \code{\link{BaseFileSection}} to be removed.}
  \item{newSection}{The \code{\link{BaseFileSection}} to replace the old one.}
  \item{append}{A \code{\link[base]{logical}}. If \code{oldSection} is not found, should
     \code{newSection} be appended anyway?}
  \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) \code{\link[base:logical]{TRUE}} if section was replaced, otherwise \code{\link[base:logical]{FALSE}}.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   For more information see \code{\link{BaseFile}}.
 }
\keyword{internal}
\keyword{methods}
