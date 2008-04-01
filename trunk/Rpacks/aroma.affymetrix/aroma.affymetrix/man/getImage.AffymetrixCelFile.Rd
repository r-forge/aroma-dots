%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% on Mon Mar 31 17:19:46 2008.
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

\usage{\method{getImage}{AffymetrixCelFile}(this, other=NULL, xrange=c(0, Inf), yrange=xrange, zrange=c(0, sqrt(2^16)), field=c("intensities", "stdvs", "pixels"), transforms=list(sqrt), interleaved=c("none", "h", "v", "auto"), ..., zoom=1, readRectFcn=NULL, verbose=FALSE)}

\arguments{
  \item{other}{An optional \code{\link{AffymetrixCelFile}} of the same chip type,
     that is used for calculating the ratio (non-logged).  Note, to get
     the log-ratios, the \code{log}() function has to be specified as
     the first transform in the \code{\link[base]{list}} of \code{transformations}.}
  \item{xrange, yrange}{\code{\link[base]{vector}}s of length two specifying the
     (x0,x1) and (y0,y1) regions to be extracted.  If \code{\link[base]{NULL}}, the
     complete regions is used.}
  \item{field}{One of the CEL file fields, i.e. \code{"intensities"},
     \code{stdvs}, or \code{pixels}.}
  \item{transforms}{A \code{\link[base]{list}} of transform \code{\link[base]{function}}s.}
  \item{interleaved}{}
  \item{...}{Additional arguments passed to
     \code{\link[aroma.affymetrix:readRawDataRectangle.AffymetrixCelFile]{*readRawDataRectangle}()}.}
  \item{zoom}{A \code{\link[base]{numeric}} scale factor in (0,+Inf) for resizing the
    imaging. If \code{1}, no resizing is done.}
  \item{readRectFcn}{A \code{\link[base]{function}} taking arguments 'xrange' and 'yrange',
    or \code{\link[base]{NULL}} for the default read function.}
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
