%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  CopyNumberSegmentationModel.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fitOne.CopyNumberSegmentationModel}
\alias{fitOne.CopyNumberSegmentationModel}
\alias{CopyNumberSegmentationModel.fitOne}
\alias{fitOne.CopyNumberSegmentationModel}
\alias{fitOne,CopyNumberSegmentationModel-method}


\title{Fits the copy-number segmentation model for one chromosome in one sample}

\description{
 Fits the copy-number segmentation model for one chromosome in one sample.

 \emph{This is an abstract method that has to be implemented in a subclass.}
}

\usage{\method{fitOne}{CopyNumberSegmentationModel}(...)}

\arguments{
  \item{data}{A \code{\link[base]{data.frame}} with columns \code{M} (log-ratio) and
     \code{x} (locus position).
  }
  \item{chromosome}{An \code{\link[base]{integer}} specifying the index of the chromosome to
     be fitted.}
  \item{...}{Additional arguments passed down to the internal fit function.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
 Returns an object returned by internal fit function.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{CopyNumberSegmentationModel}}.
}
\keyword{internal}
\keyword{methods}
