#########################################################################/**
# @RdocClass NotFilter
#
# @title "The NotFilter class"
#
# \description{
#  @classhierarchy
#
#   A NotFilter takes a set of indices from a single input and returns
#   the complementary set of indices.
# }
#
# @synopsis
#
# \arguments{
#   \item{input}{The input @see "Filter" to be connected to.}
#   \item{...}{Any arguments accepted by the @see "Filter" constructor.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
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
#    fM  <- MFilter(tma, top=0.05, col="red")
#    fT  <- TFilter(tma, top=0.05, col="blue")
#    fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
#    fSE <- NotFilter(fNotSE, visible=FALSE)
#    myFilter <- AndFilter(fM, fT, fSE, col="purple")
#
#    plot(tma, "TvsSE"); 
#    highlight(myFilter, recursive=TRUE);
# }
#
# \seealso{
#   See also the @see "ParallelFilter" class.
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("NotFilter", function(input, ...) {
  if (missing(input)) input <- NULL;
  extend(SerialFilter(input=input, ...), "NotFilter"
  )
})


setMethodS3("getIndex", "NotFilter", function(this) {
  excl <- getIndex(this$.input);
  incl <- 1:attr(excl, "max");
  idx <- setdiff(incl, excl);
  attr(idx, "max") <- attr(excl, "max");
  idx;
})



############################################################################
# HISTORY:
# 2002-02-26
# * Updated code to make use of setMethodS3's.
# 2001-07-13
# * Rearrange the class structure to contain Serial- and ParallelFilters.
# * Added some Rdoc comments.
# * Extended. Added plot paramters, recursive plotting etc.
# 2001-07-12
# * Created! Eventually I would like to have a SpeedGroupFilter, 
#   EisenFilter etc.
############################################################################
