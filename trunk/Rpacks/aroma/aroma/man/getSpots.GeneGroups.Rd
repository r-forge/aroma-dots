%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GeneGroups.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getSpots.GeneGroups}
\alias{getSpots.GeneGroups}
\alias{GeneGroups.getSpots}
\alias{getSpots.GeneGroups}
\alias{getSpots,GeneGroups-method}


\title{Gets a list of items each containing spot indices for a group}

\usage{\method{getSpots}{GeneGroups}(this, groups=NULL, unlist=FALSE, ...)}

\description{
 Gets a list of items each containing spot indices for a group.
}

\arguments{
  \item{groups}{An optional \code{\link[base]{vector}} of group indices specifying for
    which groups the spot indicies should be returned.
    If \code{\link[base]{NULL}}, all groups are used.}
}

\value{
  Returns a \code{\link[base]{list}} of length \code{nbrOfGroups(this)}.
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\examples{
SMA$loadData("mouse.data")

layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
replicates <- GeneGroups(layout)
raw <- RawData(mouse.data, layout=layout)

# Get the signal (here by default non-background corrected)
ma <- getSignal(raw)
normalizeWithinSlide(ma, "m")
normalizeAcrossSlides(ma)

print(replicates)

cat("Calculating the standardized residuals...")
idx <- getSpots(replicates)
M <- ma$M     # To improve speed!
residuals <- lapply(idx, FUN=function(i) {
  x <- M[i,]; ok <- !is.na(x)
  c <- median(x[ok]); s <- mad(x[ok], center=c)
  x <- (x-c)/s
  attr(x, "center") <- c; attr(x, "scale") <- s
  x
})
for (k in seq(idx)) M[idx[[k]],] <- residuals[[k]]
cat("done\n")

ma.res <- clone(ma)
ma.res$M <- M
subplots(ma)
for (s in seq(ma))
  plotSpatial(ma.res, slide=s)


}

\seealso{
  \code{\link[aroma:indexOf.GeneGroups]{*indexOf}()}, \code{\link[aroma:setId.GeneGroups]{*setId}()}.
  For more information see \code{\link{GeneGroups}}.
}
\keyword{internal}
\keyword{methods}
