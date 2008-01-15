%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Layout.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getPrintorder.Layout}
\alias{getPrintorder.Layout}
\alias{Layout.getPrintorder}
\alias{getPrintorder.Layout}
\alias{getPrintorder,Layout-method}


\title{Gets the order of the spots in which they were printed}

\usage{\method{getPrintorder}{Layout}(this, value=1:nbrOfSpots(this), ...)}

\arguments{
  \item{value}{The elements to be placed in the resulting \code{\link[base]{matrix}}.}
}

\description{
  Gets the order of the spots in which they were printed. The spots in column one were printed first, then the
  spots in column two and so on.
  By default the spot indices are returned.
}

\value{
  Returns a \code{\link[base]{matrix}} with spot (values) printed at the same time in the
  same column. The first spots printed are in column one and the last
  ones printed in the last column.
  Often this means that there are \code{nbrOfGrids(layout)} rows and
  \code{gridSize(layout)} columns in the matrix. However, if the
  slide was printed in say two halfs (first half of the grids are
  printed and then the second), then this is not true.
}

\section{Print order}{
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
}

\section{Different directions}{
  The most common print-order directions are \code{"row-by-row"} and
  \code{"column-by-column"}.
  In both cases, when printing a slide at each print step
  \code{nbrOfGrids(layout)} spots are printed at the same time.
  The arrayer start of spotting the first spot in \emph{each} of the
  grids. Then it cleans the print-tip heads, dries them, and go back
  to the trays to get a \emph{new} set of cDNA and prints the second
  spot in each of the grids. The second spot is to the right to
  (\code{"row-by-row"}) or below (\code{"column-by-column"}) the
  first spot. When the array gets to the end of a row (column) it
  moves on to print the next row (column) and so on until all in
  all grids have been printed.
}

\section{Print-order effects}{
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

  # No printorder specified - assumes de facto standard "row-by-row"
  print(getPrintorder(layout))

  #      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
  # [1,]    1    2    3    4    5    6    7    8    9
  # [2,]   10   11   12   13   14   15   16   17   18
  # [3,]   19   20   21   22   23   24   25   26   27
  # [4,]   28   29   30   31   32   33   34   35   36

  # Spots (1,10,19,28) were printed first, then (2,11,20,29), ...

  setPrintorder(layout, "column-by-column")
  print(getPrintorder(layout))

  #       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9]
  # [1,]    1    4    7    2    5    8    3    6    9
  # [2,]   10   13   16   11   14   17   12   15   18
  # [3,]   19   22   25   20   23   26   21   24   27
  # [4,]   28   31   34   29   32   35   30   33   36

  # Spots (1,10,19,28) were printed first, then (4,13,22,31) below, ...
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{Layout}}.
}
\keyword{internal}
\keyword{methods}
