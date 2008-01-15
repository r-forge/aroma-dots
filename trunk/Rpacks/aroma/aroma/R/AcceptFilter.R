#########################################################################/**
# @RdocClass AcceptFilter
#
# @title "The AcceptFilter class"
#
# \description{
#  @classhierarchy
#
#   An AcceptFilter asks its input for indices and let only those indices
#   through that is specified by its field \code{accept}.
# }
#
# @synopsis
#
# \arguments{
#   \item{input}{The input @see "MicroarrayData" object.}
#   \item{accept}{Vector of spot indices to be accepted.}
#   \item{...}{Other arguments accepted by the constructor of the 
#    @see "SerialFilter" class.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
#
# @author
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#  
#   ma <- getSignal(raw)
#   normalizeWithinSlide(ma, "s")
#   normalizeAcrossSlides(ma)
#  
#   tma <- as.TMAData(ma)
#  
#   # Look at the top 5\% extreme M values
#   fM  <- MFilter(tma, top=0.05, col="red")
#  
#   # And among those only look at the first 1000 spots.
#   myFilter <- AcceptFilter(fM, accept=1:1000)
#  
#   plot(tma)
#   highlight(myFilter, recursive=TRUE)
# }
#
# \seealso{
#   See also the @see "ParallelFilter" class.
# }
#*/#########################################################################
setConstructorS3("AcceptFilter", function(input, accept=NULL, ...) {
  if (missing(input)) input <- NULL;
  extend(SerialFilter(input=input, ...), "AcceptFilter", 
    accept = unlist(accept)
  )
})




#########################################################################/**
# @RdocMethod getIndex
#
# @title "Gets indices accepted by this filter"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns a @vector of @integers.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#*/######################################################################### 
setMethodS3("getIndex", "AcceptFilter", function(this) {
  incl <- getIndex(this$.input);
  idx <- intersect(incl, this$accept);
  attr(idx, "max") <- attr(incl, "max");
  idx;
})



############################################################################
# HISTORY:
# 2003-04-21
# o Added Rdoc comments.
# 2002-01-24
# * Rewritten for setClassS3 and setMethodS3.
# 2001-07-18
# * Created!
############################################################################
