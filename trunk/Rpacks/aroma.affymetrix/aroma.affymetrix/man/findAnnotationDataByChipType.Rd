%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findAnnotationDataByChipType.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{findAnnotationDataByChipType}
\alias{findAnnotationDataByChipType.default}
\alias{findAnnotationDataByChipType}

 \title{Locates an annotation data file by its chip type}

 \description{
  Locates an annotation data file by its chip type.
 }

 \usage{\method{findAnnotationDataByChipType}{default}(chipType, pattern=chipType, ..., paths=NULL, verbose=FALSE)}

 \arguments{
   \item{chipType}{A \code{\link[base]{list}} of \code{\link{AffymetrixCelSet}}.}
   \item{pattern}{A filename pattern to search for.}
   \item{...}{Arguments passed to \code{\link[affxparser]{findFiles}}.}
   \item{paths}{A \code{\link[base]{character}} \code{\link[base]{vector}} of paths to search.
     If \code{\link[base]{NULL}}, default paths are used.}
   \item{verbose}{A \code{\link[base]{logical}} or \code{\link[R.utils]{Verbose}}.}
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
