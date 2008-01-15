#########################################################################/**
# @RdocClass RejectFilter
#
# @title "The RejectFilter class"
#
# \description{
#  @classhierarchy
#
#   An RejectFilter asks its input for indices and let only those indices
#   through that is specified by its field \code{reject}.
# }
#
# @synopsis
#
# \arguments{
#   \item{input}{The input @see "Filter" to be connected to.}
#   \item{reject}{The indices to be rejected by this filter.}
#   \item{...}{Any arguments accepted by the @see "SerialFilter" constructor.}
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
#    SMA$loadData("mouse.data")
#    layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#    raw <- RawData(mouse.data, layout=layout)
#
#    ma <- getSignal(raw)
#    normalizeWithinSlide(ma, "s")
#    normalizeAcrossSlides(ma)
#
#    tma <- as.TMAData(ma)
#
#    # Look at the top 5\% extreme M values
#    fM  <- MFilter(tma, top=0.05, col="red")
#
#    # However, the first 1000 spots are not of interest
#    myFilter <- RejectFilter(fM, reject=1:1000)
#
#    plot(tma); 
#    highlight(myFilter, recursive=TRUE);
# }
#
# \seealso{
#   See also the @see "ParallelFilter" class.
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("RejectFilter", function(input, reject=NULL, ...) {
  if (missing(input)) input <- NULL;
  extend(SerialFilter(input=input, ...), "RejectFilter", 
    reject = unlist(reject)
  )
})


setMethodS3("getIndex", "RejectFilter", function(this) {
  incl <- getIndex(this$.input);
  idx <- setdiff(incl, this$reject);
  attr(idx, "max") <- attr(incl, "max");
  idx;
})



############################################################################
# HISTORY:
# 2002-02-26
# * Modified code to make use of setMethodS3's.
# 2001-07-18
# * Created!
############################################################################
