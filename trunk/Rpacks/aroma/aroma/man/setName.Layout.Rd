%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{setName.Layout}
\alias{setName.Layout}
\alias{Layout.setName}
\alias{setName.Layout}
\alias{setName,Layout-method}


\title{Sets the name of one or more spots}

\usage{\method{setName}{Layout}(this, name, index=NULL, ...)}

\arguments{
  \item{name}{A \code{\link[base]{vector}} of names.}
  \item{index}{A \code{\link[base]{vector}} of indices indicating which names to set. If
   \code{\link[base]{NULL}}, all names are set.}
}

\description{
 Sets the name of one or more spots given their indices.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.setup", "mouse.gnames"))
  layout <- as.Layout(mouse.setup, name=mouse.gnames)

  setName(layout, c("2412r", "2414r"), c(2412, 2414))

  # Get the name of spot # 2453
  getName(layout, c(2453, 2412:2417))
}

\seealso{
  \code{\link[aroma:getName.Layout]{*getName}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
