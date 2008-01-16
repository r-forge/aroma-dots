%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MicroarrayData.R
% 
% on Tue Jan 15 18:34:03 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getGeneSlideReplicateIndex.MicroarrayData}
\alias{getGeneSlideReplicateIndex.MicroarrayData}
\alias{MicroarrayData.getGeneSlideReplicateIndex}
\alias{getGeneSlideReplicateIndex.MicroarrayData}
\alias{getGeneSlideReplicateIndex,MicroarrayData-method}


\title{Returns the index of the spots in (gene,slide,replicate) indices}

\usage{\method{getGeneSlideReplicateIndex}{MicroarrayData}(this, ...)}

\description{
  Returns the index of the spots in (spot,slide) indices as
  (gene,slide,replicate) indices.
}

\value{
  Returns a \code{\link{GSRArray}} of three dimensions, gene, slide and replicate.
}

\examples{
  # Create a raw data object from the preexisting example data in
  # the sma package.
  SMA$loadData("mouse.data")

  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  # Spot 1 & 2, 3 & 4 and so on are replicates of the same gene.
  setReplicates(layout, "neighboring-pairs")

  raw <- RawData(mouse.data, layout=layout)

  # Get the signal (here by default non-background corrected)
  ma <- getSignal(raw)

  gsr <- getGeneSlideReplicateIndex(ma)

  #  Replate 1:                Replicate 2:
  #  ---------------------     ---------------------
  #  , , 1                     , , 2             
  #                                              
  #        [,1] [,2] [,3]            [,1] [,2] [,3]
  #   [1,]    1   37   73       [1,]    2   38   74
  #   [2,]    3   39   75       [2,]    4   40   76
  #   [3,]    5   41   77       [3,]    6   42   78
  #   [4,]    7   43   79       [4,]    8   44   80
  #   [5,]    9   45   81       [5,]   10   46   82
  #   [6,]   11   47   83       [6,]   12   48   84
  #   [7,]   13   49   85       [7,]   14   50   86
  #   [8,]   15   51   87       [8,]   16   52   88
  #   [9,]   17   53   89       [9,]   18   54   90
  #  [10,]   19   55   91      [10,]   20   56   92
  #  [11,]   21   57   93      [11,]   22   58   94
  #  [12,]   23   59   95      [12,]   24   60   96
  #  [13,]   25   61   97      [13,]   26   62   98
  #  [14,]   27   63   99      [14,]   28   64  100
  #  [15,]   29   65  101      [15,]   30   66  102
  #  [16,]   31   67  103      [16,]   32   68  104
  #  [17,]   33   69  105      [17,]   34   70  106
  #  [18,]   35   71  107      [18,]   36   72  108
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

\seealso{
  For more information see \code{\link{MicroarrayData}}.
}
\keyword{internal}
\keyword{methods}
