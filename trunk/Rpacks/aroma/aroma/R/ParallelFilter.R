#########################################################################/**
# @RdocClass ParallelFilter
#
# @title "The ParallelFilter class"
#
# \description{
#  @classhierarchy
#
#   A ParallelFilter is a filter that passes through indices from several
#   inputs given some criteria.
#   Examples of parallel filters are the \link{AndFilter} and the
#   @see "OrFilter", which provides the logical operators AND and OR, 
#   respectively.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{The input @see "Filter"s to be connected to.}
#   \item{cex}{The scale factor of symbols that this filter highlights.}
#   \item{col}{The color of symbols that this filter highlights.}
#   \item{pch}{The plot symbols that this filter highlights.}
#   \item{visible}{If @TRUE, the data points filtered out by this filter
#     will be highlighted, otherwise not.}
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab filters \tab The input filters. \cr
#  }
#
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
#   See also the @see "SerialFilter" class.
#   For logical filters see @see "AndFilter", @see "OrFilter", and 
#   @see "NotFilter".
#   For data filters see @see "MFilter", @see "AFilter", @see "TFilter" and 
#   @see "SEFilter".
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("ParallelFilter", function(..., cex=NULL, col=NULL, pch=NULL, visible=TRUE) {
  filters <- list(...);
  for (filter in filters) {
    if (!is.null(filter) && !inherits(filter, "Filter"))
      throw("Each argument to the constructor ParallelFilter, must be of class Filter.");
  }

  extend(Filter(cex=cex, col=col, pch=pch, visible=visible), "ParallelFilter",
    filters = filters
  )
}, abstract=TRUE)




#########################################################################/**
# @RdocMethod getInput
#
# @title "Gets the input filters"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns a @list of @see "Filter" objects.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#*/######################################################################### 
setMethodS3("getInput", "ParallelFilter", function(this) {
  this$filters;
})




#########################################################################/**
# @RdocMethod changeInput
#
# @title "Replaces input filters with other ones"
#
# @synopsis
#
# \arguments{
#  \item{recursive}{If @TRUE, the same method with the arguments are called 
#    on all input filters.}
# }
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns nothing.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#*/######################################################################### 
setMethodS3("changeInput", "ParallelFilter", function(this, newInput, oldInput=NULL, recursive=TRUE) {
  filters <- this$filters;  
  for (k in 1:length(filters)) {
    if (is.null(oldInput) || filters[[k]] == oldInput)
      filters[[k]] <- newInput
    if (recursive)
      changeInput(filters[[k]], newInput=newInput, oldInput=oldInput, recursive=TRUE);
  }
  this$filters <- filters;
})


############################################################################
# HISTORY:
# 2002-02-26
# o Modified code to make use of setMethodS3's.
# 2001-09-28
# o Changed the order of the newInput and oldInput arg's in changeInput().
# 2001-07-13
# o Rearrange the class structure to contain Serial- and ParallelFilters.
# o Added some Rdoc comments.
# o Extended. Added plot paramters, recursive plotting etc.
# 2001-07-12
# o Created! Eventually I would like to have a SpeedGroupFilter, 
#   EisenFilter etc.
############################################################################
