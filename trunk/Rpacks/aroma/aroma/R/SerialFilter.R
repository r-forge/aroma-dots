#########################################################################/**
# @RdocClass SerialFilter
#
# @title "The SerialFilter class"
#
# \description{
#  @classhierarchy
#
#   A SerialFilter is a filter that passes through indices from a single
#   input given some criteria.
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
#  \bold{Fields}
#  \tabular{rll}{
#    \tab input \tab The input. \cr
#  }
#
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
#    fM  <- MFilter(tma, top=0.05, col="red")
#    fT  <- TFilter(tma, top=0.05, col="blue")
#    fNotSE <- SEFilter(tma, range=c(-Inf,0.02), col="yellow")
#    fSE <- NotFilter(fNotSE, visible=FALSE)
#    myFilter <- AndFilter(fM, fT, fSE, col="purple")
#
#    plot(tma, "TvsSE")
#    highlight(myFilter, recursive=TRUE)
# }
#
# \seealso{
#   See also the @see "ParallelFilter" class.
# }
#*/#########################################################################
setConstructorS3("SerialFilter", function(input=NULL, ...) {
  extend(Filter(...), "SerialFilter", 
    .input = input
  )
}, abstract=TRUE)



#########################################################################/**
# @RdocMethod getInput
#
# @title "Gets all the input objects connected to the filter"
#
# @synopsis
#
# \description{
#   @get "title".
# }
#
# @author
#
# \seealso{
#   To change one or more input objects see @see "Filter.changeInput".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getInput", "SerialFilter", function(this) {
  list(this$.input);
})





#########################################################################/**
# @RdocMethod changeInput
#
# @title "Change input(s) on this filter and optionally all filters down the stream"
#
# @synopsis
#
# \description{
#   Changes some input(s) to new input(s). A filter has one or more 
#   input connections that each can be connected to for instance  another
#   filter's output or to some data objects. After having designed a
#   network of filters connected to some data objects it is sometime 
#   desireble to connect the same network of filters to another set of
#   data objects. This can be done easily by calling this method on the
#   very last filter and tell it to recursively update all other filters
#   accordingly.

# }
#
# \arguments{
#   \item{newInput}{All connections to this input will be
#     disconnected and replaced by the input object specified by 
#     this argument.}
#   \item{oldInput}{All connections matching to this input will be replaced.
#     If @NULL all connections will be replaced. Default value is
#     @NULL.}
#   \item{recursive}{If @TRUE this filter and all filters up the
#     stream will update its connections according to \code{oldInput} and
#    \code{newInput}. If @FALSE, only this filter will update its
#    input connection(s).}
# }
#
# @author
#
# \seealso{
#   To get the current set of connected inputs see @see "Filter.getInput".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("changeInput", "SerialFilter", function(this, newInput, oldInput=NULL, recursive=FALSE) {
  if (is.null(oldInput) || identical(this$.input, oldInput)) {
    this$.input <- newInput;
    if (recursive) {
      changeInput(this$.input, newInput=newInput, oldInput=oldInput, 
                                                      recursive=recursive);
    }
  }
  invisible(this);
})




############################################################################
# HISTORY:
# 2003-11-09
# o BUG FIX: Ooops, there were a really old putObject() left in the code.
# 2002-05-06
# o The SerialFilter was by mistake set to be abstract.
# 2002-02-26
# o Modified code to make use of setMethodS3's.
# 2001-09-28
# o Changed the order of the newInput and oldInput arg's in changeInput().
# 2001-07-18
# o Updated changeInput()
# 2001-07-13
# o Rearrange the class structure to contain Serial- and ParallelFilters.
# o Added some Rdoc comments.
# o Extended. Added plot paramters, recursive plotting etc.
# 2001-07-12
# o Created! Eventually I would like to have a SpeedGroupFilter, 
#   EisenFilter etc.
############################################################################
