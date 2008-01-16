%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{indexOf.Layout}
\alias{indexOf.Layout}
\alias{Layout.indexOf}
\alias{indexOf.Layout}
\alias{indexOf,Layout-method}


\title{Gets the index of one or more spots from their name or id}

\usage{\method{indexOf}{Layout}(this, name=NULL, id=name, plate=NULL, ignoreCase=TRUE, regexpr=FALSE, ...)}

\description{
 Gets the index of one or more spots from their name or id.
}

\arguments{
  \item{name}{A string \code{\link[base]{vector}} (or a single regular expression) to match
     against the names.}
  \item{id}{A string \code{\link[base]{vector}} (or a single regular expression) to match
     against the ids.}
  \item{plate}{A string \code{\link[base]{vector}} (or a single regular expression) to match
     against the plates. If it is \code{\link[base]{numeric}} it will be matched against the
     plate number as given by \code{getPlateNumber()}.}
  \item{ignoreCase}{If \code{\link[base:logical]{TRUE}}, the matching is not case sensitive,
     otherwise it is.}
  \item{regexpr}{If \code{\link[base:logical]{TRUE}}, regular expression matching is used,
     otherwise plain string comparison is used.}
  At least one of the arguments \code{name} and \code{id} must be given.
  If both are given, indices of spots that match \emph{either} the
  \code{name} search pattern \emph{or} the \code{id} search pattern.
}

\value{
  Returns the \code{\link[base]{vector}} of indices of matched names or ids. Returns
  \code{numeric(0)} if no matches were found.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
  SMA$loadData(c("mouse.setup", "mouse.gnames"))
  layout <- as.Layout(mouse.setup, id=mouse.gnames)

  # Get the index of spots with id "54" and "232".
  indexOf(layout, id=c("54", "232"))
  # [1] 54 232

  # Get the index of all spots with id beginning with "120" and
  # having at least four characters.
  indexOf(layout, id="^120.+", regexpr=TRUE)
  # [1] 1200 1201 1202 1203 1204 1205 1206 1207 1208 1209
}

\seealso{
  For more help on regular expressions see
  \code{\link[base]{grep}}() and \code{\link[base]{apropos}}.
  \code{\link[aroma:getName.Layout]{*getName}()}.
  \code{\link[aroma:getId.Layout]{*getId}()}.
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
