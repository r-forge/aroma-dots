\name{findUnit}
\alias{findUnit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ }
\description{Find unit name for a given group name.}
\usage{
findUnit(cdf, names, groupList = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{ cdf to search }
  \item{names}{ vector of group names to find }
}

\value{
    vector of unit numbers (cdf dependent) corresponding to each group name. If there was no group name assigned element NA; if multiple matches, gives the first match with warning.
}
\author{Elizabeth Purdom}
