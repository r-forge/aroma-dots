%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.PLOT.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{writeImage.AffymetrixCelFile}
\alias{writeImage.AffymetrixCelFile}
\alias{AffymetrixCelFile.writeImage}
\alias{writeImage.AffymetrixCelFile}
\alias{writeImage,AffymetrixCelFile-method}


\title{Writes a spatial image of the signals in the CEL file}

\description{
 Writes a spatial image of the signals in the CEL file.
}

\usage{\method{writeImage}{AffymetrixCelFile}(this, filename=NULL, fullname=getFullName(this), tags=c("*", "gray"), imgFormat="png", path=NULL, ..., verbose=FALSE)}

\arguments{
  \item{...}{Arguments passed to \code{\link[aroma.affymetrix:getImage.AffymetrixCelFile]{*getImage}()}.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns (invisibly) the Image written.
}

\examples{\dontrun{
  df <- getFile(ds, 1)
  writeImage(df, tags="gray", palette=gray.colors(256), xrange=c(0,200))
  writeImage(df, tags="heat", palette=heat.colors(256), xrange=c(0,200))
  writeImage(df, tags="terrain", palette=terrain.colors(256), xrange=c(0,200))
  writeImage(df, tags="topo", palette=topo.colors(256), xrange=c(0,200))
  writeImage(df, tags="cm", palette=cm.colors(256), xrange=c(0,200))
  writeImage(df, tags="rainbow", palette=rainbow(256), xrange=c(0,200))
}}

\author{Ken Simpson (ksimpson[at]wehi.edu.au).}

\seealso{
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
