%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getLabel.MicroarrayData}
\alias{getLabel.MicroarrayData}
\alias{MicroarrayData.getLabel}
\alias{getLabel.MicroarrayData}
\alias{getLabel,MicroarrayData-method}


\title{Gets the label of one field}

\usage{\method{getLabel}{MicroarrayData}(this, field, ...)}

\arguments{
 \item{field}{A field name to get the label for.}
}

\description{
  Gets the label of one field.
  Field labels are for instance used as default value is many plot
  functions.
  The labels can be set with the \code{\link[aroma:setLabel.MicroarrayData]{*setLabel}()} method.
}

\value{
  Returns the label. If the field was not found, the field is returned.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  print(getLabel(raw, "R"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
 \code{\link[aroma:setLabel.MicroarrayData]{*setLabel}()}.
 See \code{\link[grDevices]{plotmath}} how expressions can be used to annotate
 graphs.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
