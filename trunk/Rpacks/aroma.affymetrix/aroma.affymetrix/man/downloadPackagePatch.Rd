%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  downloadPackagePatch.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{downloadPackagePatch}
\alias{downloadPackagePatch.default}
\alias{downloadPackagePatch}

 \title{Download a package patch}

 \description{
  Download a package patch from the package online reprocitory.
 }

 \usage{\method{downloadPackagePatch}{default}(pkgName, version=NULL, url=NULL, verbose=FALSE, ...)}

 \arguments{
   \item{pkgName}{The name of the package to be patched."}
   \item{version}{A \code{\link[base]{character}} string.}
   \item{url}{The root URL from where to download the patch.}
   \item{...}{Not used.}
 }

 \value{
   Returns (invisibly) the pathnames of the files downloaded.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}


\keyword{internal}
