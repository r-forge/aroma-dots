%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% on Tue Jan 15 18:32:39 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{toPrintorderMatrix.Layout}
\alias{toPrintorderMatrix.Layout}
\alias{Layout.toPrintorderMatrix}
\alias{toPrintorderMatrix.Layout}
\alias{toPrintorderMatrix,Layout-method}


\title{Gets a matrix of spot indices in the order they were printed}

\usage{\method{toPrintorderMatrix}{Layout}(this, value=seq(nbrOfSpots(this)), ...)}

\arguments{
  \item{value}{The elements to be placed in the resulting matrix.}
}

\description{
  Gets a matrix of spot indices in the order they were printed. The
  spots in column one were printed first, then the spots in column two
  and so on. By default the spot indices are returned.
}

\value{
  Returns a \code{\link[base]{matrix}} with \code{nbrOfGrids(layout)} rows and
  \code{gridSize(layout)} columns.
}

\details{
  When printing a slide at each print step \code{nbrOfGrids(layout)}
  spots are printed at the same time. The arrayer start of spotting the
  first spot in \emph{each} of the grids. The it cleans the print-tip
  heads, dries them, and go back to the trays to get a \emph{new} set of
  cDNA and prints the second spot in each of the grids. The second spot
  is to the right of the first spot. When the array gets to the end of
  a row it has printed the first line of each grid. The arrayer now moves
  on to print the spots in line two, three and so on until all lines in
  all grids have been printed.\cr

  The printing of a microarray is time consuming and often several
  microarray slides are printed at the same time, since it is even more
  time consuming to switch between the trays. When printing several
  microarrays at the same time, the arrayer prints the first spot in
  all grids on \emph{all} slides, before moving on to the second spot.
  For a example, printing a batch of 100 slides with 6384 spots in 4x4
  grids takes about 15 hours to print including manual work to switch
  trays etc. Each grid contains 19*21 spots, i.e. the arrayer has to put
  down the print tips 399 times on each slide, and in total 39900 times.
  This is about 44 put-downs a minute. It takes about 45-50 minutes
  to finish one row of spots.

  An important factor for the quality of the printed spots is the
  temperature and the humidity. Too high temperature and humidity
  tends to produce too large spots that can even overlap [1]. If there
  is no automatic control for temperature and humidity, the quality of
  the spots could vary a lot between the spots printed during a 15 hours
  printing process. With a varying printing climate we should expect to
  see a variating of the quality of the spots along the order of which
  the spots are printed.
  The variation of temperature and humidity varies approximately in the
  time scale of hours. As it takes about 45-50 minutes to print a row
  of spots, we should therefore expect to see such a variation between
  the rows in the grids.
}

\references{
  [1] Microarrays in Three Easy Steps, Priti Hedge, The Institute for
      Genomic Research, 200?.
}

\examples{
  layout <- Layout(2,2, 3,3)

  print(toPrintorderMatrix(layout))

  #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
  # [1,]    1    2    3    4    5    6    7    8    9
  # [2,]   10   11   12   13   14   15   16   17   18
  # [3,]   19   20   21   22   23   24   25   26   27
  # [4,]   28   29   30   31   32   33   34   35   36

  # Spot 1, 10, 19 and 28 were printed first, then 2, 11, 20, and 29...
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
