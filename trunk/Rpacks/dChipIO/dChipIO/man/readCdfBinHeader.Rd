%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readCdfBinHeader.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readCdfBinHeader}
\alias{readCdfBinHeader}


\title{Reads the file header of a dChip CDF.bin file}

\description{
  Reads the file header of a dChip CDF.bin file.
}

\usage{readCdfBinHeader(con, ...)}

\arguments{
  \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} filename.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} structure containing the file header.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To read the CDF.bin file data, see \code{\link{readCdfBin}}().
}



\keyword{file}
\keyword{IO}
