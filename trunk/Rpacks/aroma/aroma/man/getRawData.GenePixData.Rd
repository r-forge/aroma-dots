%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenePixData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getRawData.GenePixData}
\alias{getRawData.GenePixData}
\alias{GenePixData.getRawData}
\alias{getRawData.GenePixData}
\alias{getRawData,GenePixData-method}


\title{Gets the raw intensites from the GPR structure}

\description{
 Extracts the red and green spot intensitites (both signal and background)
 from the GenePixData object and returns a \code{\link{RawData}} object.

 \emph{Note: From com.braju.sma v0.59, this method returns the
       *median* foreground and background estimates. Previous
       versions returned the *mean* of dito.}
}

\usage{\method{getRawData}{GenePixData}(this, slides=NULL, fg="median", bg="median", channels=1:2, ...)}

\arguments{
  \item{slides}{Specifying which slides to be extracted. If \code{\link[base]{NULL}}, all slides are considered.}
  \item{fg}{If \code{"mean"}, the mean foreground intensities are returned. If \code{"median"}, the median foreground intensities are returned.}
  \item{bg}{If \code{"mean"}, the mean background intensities are returned. If \code{"median"}, the median background intensities are returned.}
  \item{channels}{A \code{\link[base]{vector}} of length two specifying which two channels
    (wavelengths) to be extracted in case the data contains more than two
     channels (colours).}
}

\value{Returns a \code{\link{RawData}} object containing the specified slides.}

\details{
  The R and Rb channels will come from the F635* and B635* fields, and
  the G and Gb channels will come from the F532* and B532* fields.
  To swap the channels just use dyeSwap().

  From com.braju.sma v0.30 the \code{dye.swap} argument has been removed.
  Any dye swapping should be done bye using the dyeSwap() method.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\references{
  GenePix File Formats,
  \url{http://www.axon.com/GN_GenePix_File_Formats.html}
}

\examples{
  gpr <- GenePixData$read("gpr123.gpr", path=system.file("data-ex", package="aroma"))

  # Gets the foreground and the background
  raw <- getRawData(gpr)
}

\seealso{
  For more information see \code{\link{GenePixData}}.
}
\keyword{internal}
\keyword{methods}
