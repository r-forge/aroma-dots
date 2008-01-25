\name{findUnit}
\alias{findUnit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Unit from Group Name }
\description{Find unit name for a given group name.}
\usage{
findUnit(cdf, names, returnNames = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cdf}{ cdf to search }
  \item{names}{ vector of group names to find }
  \item{returnNames}{ logical, whether return unit names as well as number (may have problems)}
}

\value{
    vector of unit numbers (cdf dependent) corresponding to each group name. If there was no group name assigned element NA; if multiple matches, gives the first match with warning.
}
\author{Elizabeth Purdom}
