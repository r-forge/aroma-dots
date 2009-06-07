%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{read.Layout}
\alias{read.Layout}
\alias{Layout.read}
\alias{read.Layout}
\alias{read,Layout-method}


\title{Reads layout information from a tab-delimited file}

\description{
 Static method that reads layout information from a tab-delimited file.
}

\synopsis{read.Layout(this, filename, path=NULL, header=TRUE, sep="	", quote="", ...)}
\usage{Layout$read(filename, path=NULL, header=TRUE, sep="	", quote="", ...)}

\arguments{
  \item{filename}{The name of the file.}
  \item{path}{Optional path where the data should be written.}
  \item{sep}{Separator \code{\link[base]{character}} between columns.}
  \item{header}{If \code{\link[base:logical]{TRUE}} column headers are written, otherwise not.}
  \item{...}{Other arguments accepted by \code{read.table()}.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  print(layout)
}

\seealso{
  For writing a Layout object to a file \code{\link[aroma:write.Layout]{*write}()}.
  See also \code{read.table()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
