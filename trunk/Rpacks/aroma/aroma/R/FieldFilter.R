#########################################################################/**
# @RdocClass FieldFilter
#
# @title "The FieldFilter class"
#
# \description{
#  @classhierarchy
#
#   A FieldFilter provides methods to extract indices from
#   @see "MicroarrayData" objects given some criteria on one of the 
#   \emph{fields}. Possible criterias on this field are \emph{top}, 
#   and \emph{range}.
#   The criteria \emph{top} filters out the top \emph{N} fraction 
#   (or number if \emph{N>1}).
#   The criteria \emph{range} filters accepts the spots with a field value
#   within the given range.
# }
#
# @synopsis
#
# \arguments{
#   \item{mad}{A @see "MicroarrayData" object to be filtered.}
#   \item{field}{The field (@character name) to be filtered.}
#   \item{top, bottom}{If specified, to top (bottom) values are filtered out. 
#    If an @integer one or greater, that number of indicies will be passed.
#    If a @numeric between zero and one, that ratio will be passed.}
#   \item{range}{The range of values to be passed.}
#   \item{absolute.values}{If @TRUE, absolute values are filtered, 
#    otherwise not.}
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
#    \tab range \tab The range criteria. \cr
#    \tab top \tab The top criteria. \cr
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
#    plot(tma, "TvsSE"); 
#    highlight(myFilter, recursive=TRUE);
# }
#
# \seealso{
#   See also the @see "Filter" class.
#   For logical filters see @see "AndFilter", @see "OrFilter", and 
#   @see "NotFilter".
#   For data filters see @see "MFilter", @see "AFilter", @see "TFilter" and 
#   @see "SEFilter".
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("FieldFilter", function(mad, field, bottom=NULL, top=NULL, range=NULL, absolute.values=TRUE, cex=NULL, col=NULL, pch=NULL, visible=TRUE) {
  ok <- !missing(mad);
  if (ok && !inherits(mad, "MicroarrayData"))
    throw("Argument mad in the constructor FieldFilter, must be of class MicroarryData.");
  if (!ok) mad <- NULL;

  extend(SerialFilter(input=mad, cex=cex, col=col, pch=pch, visible=visible), "FieldFilter", 
    field = if (ok) field else NULL,
    range = range,
    bottom = bottom,
    top = top,
    absolute.values = absolute.values
  )
})


#########################################################################/**
# @RdocMethod getIndex
#
# \title{Gets indices accepted by this filter}
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
setMethodS3("getIndex", "FieldFilter", function(this) {
  field <- this$field;
  mad <- this$.input;
  x <- extract(mad, field=field);
  getIndex0.FieldFilter(this, x);
})


setMethodS3("getIndex0", "FieldFilter", function(this, x) {
  require(R.basic) || throw("Package R.basic was not found!"); # top()

  mad <- this$.input;
  
  x <- as.matrix(x);  # Has to be done since 'x' can be a data.frame.
  if (this$absolute.values)
    x <- abs(x);
  
  if (!is.null(this$bottom))
    idx <- which(top(-x, n=this$top))
  else if (!is.null(this$top))
    idx <- which(top(x, n=this$top))
  else
  idx <- seq(x);
    
  range <- this$range;
  if (!is.null(range))
    idx <- intersect(idx, which(x >= range[1] & x <= range[2]));
  attr(idx, "max") <- size(mad); # Max index.
  idx;
}, protected=TRUE);





############################################################################
# HISTORY:
# 2003-05-04
# o Updated the Rdoc with arguments.
# 2003-04-21
# o Added Rdocs.
# 2002-02-26
# * Updated code to make use of setMethodS3's.
# 2002-01-24
# * Renamed all get() to extract()
# 2001-07-15
# * Removed getIndex.SEFilter since SE now is a true field of TMAData.
# 2001-07-13
# * Rearrange the class structure to contain Serial- and ParallelFilters.
# * Added some Rdoc comments.
# * Extended. Added plot paramters, recursive plotting etc.
# 2001-07-12
# * Created! Eventually I would like to have a SpeedGroupFilter, 
#   EisenFilter etc.
############################################################################
