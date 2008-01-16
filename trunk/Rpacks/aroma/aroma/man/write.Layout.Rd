%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{write.Layout}
\alias{write.Layout}
\alias{Layout.write}
\alias{write.Layout}
\alias{write,Layout-method}


\title{Writes the layout information to a tab-delimited file}

\description{
 Writes the layout information to a tab-delimited file.
}

\usage{\method{write}{Layout}(this, filename, path=NULL, overwrite=FALSE, row.names=FALSE, sep="	", quote=FALSE, ...)}

\arguments{
  \item{filename}{The name of the file.}
  \item{path}{Optional path where the data should be written.}
  \item{row.names}{If \code{\link[base:logical]{TRUE}} row names are written, otherwise not.}
  \item{sep}{Separator \code{\link[base]{character}} between columns.}
  \item{...}{Other arguments accepted by \code{write.table()}.}
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  layout <- Layout(3,2, 2,3)
  tmpfile <- tempfile()
  write(layout, tmpfile)
  file.show(tmpfile)
  unlink(tmpfile)
}

\seealso{
  For read a Layout object from file see \code{\link[aroma:read.Layout]{*read}()}.
  See also \code{write.table()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
