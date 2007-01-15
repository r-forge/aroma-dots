%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
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

\usage{\method{getImage}{AffymetrixCelFile}(this, xrange=c(0, Inf), yrange=c(0, Inf), takeLog=TRUE, interleaved=FALSE, lim=c(-Inf, Inf), field=c("intensities", "stdvs", "pixels"), col=NULL, outlierCol="white", main=getName(this), ...)}

\arguments{
  \item{}{}
}

\value{
  Returns an Image object.
}

\section{Details}{
}

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
