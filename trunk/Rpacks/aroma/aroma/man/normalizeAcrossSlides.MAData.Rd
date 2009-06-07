%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.NORM.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeAcrossSlides.MAData}
\alias{normalizeAcrossSlides.MAData}
\alias{MAData.normalizeAcrossSlides}
\alias{normalizeAcrossSlides.MAData}
\alias{normalizeAcrossSlides,MAData-method}


\title{Normalizes across slides}

\usage{\method{normalizeAcrossSlides}{MAData}(this, slides=NULL, newMAD=NULL, ...)}

\arguments{
  \item{slides}{The set slides to be used and to be normalized. If
    \code{\link[base]{NULL}} all slides will be normalized.}
  \item{newMAD}{After the normalization all slides will have an maximum
    absolute deviation of \code{newMAD}. If \code{\link[base]{NULL}}, the result will
    be the same as if \code{newMAD} would be set to the geometrical mean
    of each individual slide's MAD.}
}

\description{
  Normalizes across some or all slides. After doing within-slide normalization
  (see \code{\link[aroma:normalizeWithinSlide.MAData]{*normalizeWithinSlide}()}, an across-slide normalization must
  be performed before the data on the different slides can be compared with
  each. Across-slide normalization scales the log ratios (M) for all slide so
  each slide gets the same log-ratio deviation based on the robust deviation
  measure Maximum Absolute Deviation (MAD).

  If one would like to normalize the deviation across \emph{different data set},
  i.e. different \code{\link{MAData}} objects, one can make use of the argument
  \code{newMAD}, which forces the slides to get a specific median absolute
  deviation value.

  Note, in the case where one set of slides comes from one type of experimental
  setup and a second set of slides comes from another setup, and they are stored
  in the same \code{MAData} object, these two groups of slides \emph{can} be
  normalized together using one (in other words, you do \emph{not} have to
  normalize the two groups seperately and the rescale them with \code{newMAD}).

  Also note that it is only the log ratios, M, are affected by the
  normalization, i.e. the log intensities, A, are \emph{not} changed.
}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  layout(matrix(1:9, ncol=3, byrow=TRUE))

  # Plot data before within-slide normalization
  for (k in 1:3)
    boxplot(ma, groupBy="printtip", slide=k, main=paste("No normalization - #", k, sep=""))

  # Plot data after scaled print-tip normalization
  normalizeWithinSlide(ma, "s")
  for (k in 1:3)
    boxplot(ma, groupBy="printtip", slide=k, main=paste("Within-slide norm. - #", k, sep=""))

  # Plot data after across-slide normalization
  normalizeAcrossSlides(ma)
  for (k in 1:3)
    boxplot(ma, groupBy="printtip", slide=k, main=paste("Across-slide norm. - #", k, sep=""))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For an detailed explanation of the robust deviation measure MAD
  (Median Absolute Deviation) see \code{\link[stats]{mad}}.
  For within-slide normalization see \code{\link[aroma:normalizeWithinSlide.MAData]{*normalizeWithinSlide}()}.
  For more information see \code{\link{MAData}}.
}
\keyword{internal}
\keyword{methods}
