%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getId.Layout}
\alias{getId.Layout}
\alias{Layout.getId}
\alias{getId.Layout}
\alias{getId,Layout-method}


\title{Gets the id of one or more spots}

\usage{\method{getId}{Layout}(this, index=NULL, ...)}

\arguments{
  \item{index}{A \code{\link[base]{vector}} of indices indicating which ids to set. If
   \code{\link[base]{NULL}}, all ids are set.}
}

\value{
  Returns the \code{\link[base]{vector}} of ids.
}

\description{
 Gets the id of one or more spots given their indices.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.setup", "mouse.gnames"))
  layout <- as.Layout(mouse.setup, id=mouse.gnames)

  # Get the id of spot # 2453 and 2412:2417.
  getId(layout, c(2453, 2412:2417))
}

\seealso{
  \code{\link[aroma:indexOf.Layout]{*indexOf}()},
  \code{\link[aroma:setId.Layout]{*setId}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
