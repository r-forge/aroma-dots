%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{extract.MicroarrayData}
\alias{extract.MicroarrayData}
\alias{MicroarrayData.extract}
\alias{extract.MicroarrayData}
\alias{extract,MicroarrayData-method}


\title{Gets a subset of data specified by fields, slides and/or spot indices}

\usage{\method{extract}{MicroarrayData}(this, fields=NULL, slides=NULL, index=NULL, ...)}

\description{
  Gets a subset of data specified by fields, slides and/or spot indices.
  Returns the specified subset of the microarray data as data frame with
  a first column specifying the slide.
  If the \code{fields} argument is not specified, all fields are considered.
  If the \code{slides} argument is not specified, all slides are considered.
  All data elements returned matches by field (name or index) AND by slide
  number AND by spot index.
}

\arguments{
  \item{fields}{The field names or indices to be returned. If \code{\link[base]{NULL}},
    all fields are included.}
  \item{slides}{The slides that should be included. If \code{\link[base]{NULL}},
    all slides are included.}
  \item{index}{The spot indices that should be included. If \code{\link[base]{NULL}},
    all spots are included.}
  Any of the above arguments are optional.
}

\value{Returns a \code{\link[base]{data.frame}}.}

\details{
  This method was formerly named \code{get()}, but since this is a
  insecure name of a method it was renamed. In addition to the fields
  as of \code{getFieldNames()}, there are also three special fields:
  1) \code{slide}, 2) \code{spot} and 3) \code{gene}. These fields are
  not accessable fields per se, but they are calculated on request.
  The field \code{gene} is only returned if there are within-slide
  replicates defined (see class \code{\link{Layout}}).
  For instance, \code{as.data.frame()} asks for these fields internally.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  # Get the log ratios (M) for spots 4-7 in the slide 2,3 and 4.
  extract(ma, c("slide", "M"), slides=2:4, index=4:7)
  # Gives:
  #      slide           M
  #   1      2 -0.41129383
  #   2      2 -0.44570800
  #   3      2 -0.16409736
  #   4      2 -0.22462037
  #   5      3 -0.51599488
  #   6      3 -0.30718292
  #   7      3 -0.54794073
  #   8      3 -0.42497673
  #   9      4  0.51019368
  #   10     4  0.31210389
  #   11     4  0.08827923
  #   12     4 -0.12370050
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:as.data.frame.MicroarrayData]{*as.data.frame}()}.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
