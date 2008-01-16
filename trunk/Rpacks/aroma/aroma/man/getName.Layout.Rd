%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getName.Layout}
\alias{getName.Layout}
\alias{Layout.getName}
\alias{getName.Layout}
\alias{getName,Layout-method}


\title{Gets the name of one or more spots}

\usage{\method{getName}{Layout}(this, index=NULL, ...)}

\arguments{
  \item{index}{A \code{\link[base]{vector}} of indices indicating which names to set. If
   \code{\link[base]{NULL}}, all names are set.}
}

\value{
  Returns the \code{\link[base]{vector}} of names.
}

\description{
 Gets the name of one or more spots given their indices.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.setup", "mouse.gnames"))
  layout <- as.Layout(mouse.setup, name=mouse.gnames)

  # Get the name of spot # 2453 and 2412:2417.
  getName(layout, c(2453, 2412:2417))
}

\seealso{
  \code{\link[aroma:indexOf.Layout]{*indexOf}()},
  \code{\link[aroma:setName.Layout]{*setName}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
