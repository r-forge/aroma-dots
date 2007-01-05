%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AffymetrixCelFile.R
% 
% on Sat Jan 06 01:47:01 2007.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{writeSpatial.AffymetrixCelFile}
\alias{writeSpatial.AffymetrixCelFile}
\alias{AffymetrixCelFile.writeSpatial}
\alias{writeSpatial.AffymetrixCelFile}
\alias{writeSpatial,AffymetrixCelFile-method}


\title{Writes a spatial image of the CEL data to an image file}

\description{
 Writes a spatial image of the CEL data to an image file.
}

\usage{\method{writeSpatial}{AffymetrixCelFile}(this, filename=sprintf(fmtstr, getName(this)), path=filePath("figures", getChipType(getCdf(this))), fmtstr="\%s-spatial.png", ..., field=c("intensities", "stdvs", "pixels"), transform=sqrt, zlim=NULL, verbose=FALSE)}

\arguments{
  \item{filename}{The filename of the image file.}
  \item{path}{The path to the image file.}
  \item{...}{Additional arguments passed to \code{\link[aroma.affymetrix:readRectangle.AffymetrixCelFile]{*readRectangle}()}.}
  \item{field}{The data field to be displayed.}
  \item{verbose}{A \code{\link[base]{logical}} or a \code{\link[R.utils]{Verbose}} object.}
}

\value{
  Returns (invisibly) the pathname of the generated image file.
}

\section{Details}{
  If any other image format than the pixmap *.pgm format is requested,
  an external image coverted is used to convert from the PGM format.
  You may specify the absolute pathname to such an image converted by:
  \code{options(imageConverter="C:/Program Files/ImageMagick/bin/convert")}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma.affymetrix:image270.AffymetrixCelFile]{*image270}()}.
  For more information see \code{\link{AffymetrixCelFile}}.
}


\keyword{internal}
\keyword{methods}
\keyword{IO}
