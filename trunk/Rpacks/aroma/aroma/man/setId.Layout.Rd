%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setId.Layout}
\alias{setId.Layout}
\alias{Layout.setId}
\alias{setId.Layout}
\alias{setId,Layout-method}


\title{Sets the id of one or more spots}

\usage{\method{setId}{Layout}(this, id, index=NULL, ...)}

\arguments{
  \item{id}{A \code{\link[base]{vector}} of ids.}
  \item{index}{A \code{\link[base]{vector}} of indices indicating which ids to set. If
   \code{\link[base]{NULL}}, all ids are set.}
}

\description{
 Sets the id of one or more spots given their indices.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.setup", "mouse.gnames"))
  layout <- as.Layout(mouse.setup, id=mouse.gnames)

  setId(layout, c("2412r", "2414r"), c(2412, 2414))

  # Get the id of spot # 2453
  getId(layout, c(2453, 2412:2417))
}

\seealso{
  \code{\link[aroma:getId.Layout]{*getId}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
