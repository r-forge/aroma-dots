%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{equals.MicroarrayData}
\alias{equals.MicroarrayData}
\alias{MicroarrayData.equals}
\alias{equals.MicroarrayData}
\alias{equals,MicroarrayData-method}


\title{Compares a MicroarrayData object with another Object}

\usage{\method{equals}{MicroarrayData}(this, other, ...)}

\description{
 Compares a MicroarrayData object with another Object and returns
 \code{\link[base:logical]{TRUE}} if they are equal. For the MicroarrayData object \code{raw}
 to be equal to another Object \code{obj} it must be true that:

 \item{1}{\code{obj} is of class MAData,}
 \item{2}{\code{size(raw) == size(obj)},}
 \item{3}{\code{equals(getLayout(raw), getLayout(obj))} is true, and}
 \item{4}{the Euclidian distance between the each field in the two objects
   is small.}
}

\examples{
  # The option 'dataset' is used to annotate plots.
  options(dataset="sma:MouseArray")

  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  ma <- getSignal(raw)

  ma2 <- clone(ma)
  print(equals(ma, ma2))   # TRUE

  layout2 <- Layout(4,4,21,19)
  setLayout(ma2, layout2)
  print(equals(ma, ma2))   # FALSE
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
