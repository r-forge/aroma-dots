%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  ScanAlyze20Data.R
% 
% on Tue Jan 15 18:36:41 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{write.ScanAlyze20Data}
\alias{write.ScanAlyze20Data}
\alias{ScanAlyze20Data.write}
\alias{write.ScanAlyze20Data}
\alias{write,ScanAlyze20Data-method}


 \title{Writes a ScanAlyze v2.0 Results Data file (NOT IMPLEMENTED)}

 \description{
  Writes a ScanAlyze v2.0 Results Data file (NOT IMPLEMENTED).

  Not implemented yet!
 }

 \usage{\method{write}{ScanAlyze20Data}(this, filename, path=NULL, slides=NULL, row.names=FALSE, ..., verbose=FALSE)}

 \arguments{
   \item{filename}{The filename of the file to be read.}
   \item{path}{The path to the file.}
   \item{slides}{An \code{\link[base]{integer}} specifying which slides to be written to file. Currently, only one slide at the time can be written.}
   \item{...}{Arguments passed to \code{write.table}.}
 }

 \value{Returns nothing.}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   To read one or more ScanAlyze v2.0 Results files see \code{\link[aroma:read.ScanAlyze20Data]{*read}()}.
   For more information see \code{\link{ScanAlyze20Data}}.
 }
\keyword{internal}
\keyword{methods}
