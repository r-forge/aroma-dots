%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{equals.Layout}
\alias{equals.Layout}
\alias{Layout.equals}
\alias{equals.Layout}
\alias{equals,Layout-method}


\title{Checks if a Layout object is equals to some other object}

\usage{\method{equals}{Layout}(this, obj, ...)}

\arguments{
 \item{obj}{The other object for which this object should be compared to.}
}

\description{
 Checks if a \code{Layout} object is equal to some other object, which
 normally is another \code{Layout} object. It could however also be a
 list with the same fields as a \code{Layout} object.
 \cr
 A layout is equal to another layout if it has 1) the same number of
 rows and columns of grids as the other object, 2) the same number of
 rows and columns of spots as the other object.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
   layout1 <- Layout(4,4, 18,18)
   layout2 <- Layout(4,4, 18,18)
   layout3 <- Layout(4,4, 16,16)
   equals(layout1, layout1)   # TRUE (of course)
   equals(layout1, layout2)   # TRUE
   equals(layout2, layout1)   # TRUE (symmetric)
   equals(layout1, layout3)   # FALSE
   equals(layout2, layout3)   # FALSE
}
\keyword{internal}
\keyword{methods}
