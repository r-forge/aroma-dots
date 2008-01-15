%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setLabel.MicroarrayData}
\alias{setLabel.MicroarrayData}
\alias{MicroarrayData.setLabel}
\alias{setLabel.MicroarrayData}
\alias{setLabel,MicroarrayData-method}


\title{Sets the label of one field}

\usage{\method{setLabel}{MicroarrayData}(this, field, value, ...)}

\arguments{
 \item{field}{A field name to which the label should be assigned.}
 \item{label}{The new label.}
}

\description{
  Sets the label of one field (only). If the field does not exists the
  label will stored anyway under a virtual field, which could be useful
  to store different labels.
}

\examples{
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  setLabel(raw, "R", "red")
  print(getLabel(raw, "R"))
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  \code{\link[aroma:getLabel.MicroarrayData]{*getLabel}()}.
  See \code{\link[grDevices]{plotmath}} how expressions can be used to annotate
  graphs.
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
