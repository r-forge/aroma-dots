%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% on Tue Nov 20 12:39:02 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getImage.AffymetrixCelFile}
\alias{getImage.AffymetrixCelFile}
\alias{AffymetrixCelFile.getImage}
\alias{getImage.AffymetrixCelFile}
\alias{getImage,AffymetrixCelFile-method}


\title{Creates an RGB Image object from a CEL file}

\description{
 Creates an RGB Image object from a CEL file.
}

\usage{\method{getImage}{AffymetrixCelFile}(this, xrange=c(0, Inf), yrange=xrange, zrange=c(0, sqrt(2^16)), field=c("intensities", "stdvs", "pixels"), transforms=list(sqrt), interleaved=c("none", "h", "v", "auto"), ..., zoom=NULL, verbose=FALSE)}

\arguments{
  \item{xrange, yrange}{}
  \item{field}{}
  \item{transforms}{A \code{\link[base]{list}} of transform functions.}
  \item{interleaved}{}
  \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:getRectangle.AffymetrixCelFile]{*getRectangle}()}.}
  \item{zoom}{}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns an Image object.
}

\section{Details}{
}

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}

\seealso{
  \code{\link[EBImage]{Image}}.
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
