%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readDcpHeader.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readDcpHeader}
\alias{readDcpHeader}


\title{Reads the file header of a dChip DCP file}

\description{
  Reads the file header of a dChip DCP file.
}

\usage{readDcpHeader(con, ...)}

\arguments{
  \item{con}{A \code{\link[base:connections]{connection}} or a \code{\link[base]{character}} filename.}
  \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{list}} structure containing the file header.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  To read also the DCP file data, see \code{\link{readDcp}}().
}



\keyword{file}
\keyword{IO}
