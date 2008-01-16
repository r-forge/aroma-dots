%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{set.Layout}
\alias{set.Layout}
\alias{Layout.set}
\alias{set.Layout}
\alias{set,Layout-method}


\title{Sets the layout}

\usage{\method{set}{Layout}(this, ngrid.r=NULL, ngrid.c=NULL, nspot.r=NULL, nspot.c=NULL, name=NULL, id=NULL, ...)}

\arguments{
  \item{nspot.r}{The number of rows of spots per grid.}
  \item{nspot.c}{The number of columns of spots per grid.}
  \item{ngrid.r}{The number of rows of grids per slide.}
  \item{ngrid.c}{The number of columns of grids per slide.}
  \item{name}{A \code{\link[base]{vector}} if names for the spot.}
  \item{id}{A \code{\link[base]{vector}} if ids for the spot.}
}

\description{
 Sets the layout by either 1) explicitly setting the number of rows and
 columns of grids and the number of rows and columns of spots within each
 grid, or 2) by giving another \code{Layout} object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   layout1 <- Layout(4,4, 18,18)
   layout2 <- Layout()
   set(layout2, 4,4, 16,16)   # Alternative 1
   set(layout2, layout1)      # Alternative 2
}
\keyword{internal}
\keyword{methods}
