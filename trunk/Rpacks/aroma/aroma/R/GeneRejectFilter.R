#########################################################################/**
# @RdocClass GeneRejectFilter
#
# @title "The GeneRejectFilter class"
#
# \description{
#  @classhierarchy
#
#   An GeneRejectFilter asks its input for indices and let only those
#   indices through that do \emph{not} corresponds to a given set of genes.
# }
#
# @synopsis
#
# \arguments{
#   \item{input}{The input @see "Filter" to be connected to.}
#   \item{layout}{A @see "Layout" object.}
#   \item{genes}{Gene}.
#   \item{ids}{}.
#   \item{names}{}.
#   \item{...}{Any arguments accepted by the @see "RejectFilter" constructor.}
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
#    # Look at the top 5\% extreme M values
#    fM  <- MFilter(tma, top=0.05, col="red")
#
#    # However, the blanks, e.g. empty spots, are not of interest
#    myFilter <- GeneRejectFilter(fM, layout=layout, ids="BLANK", col="orange")
#
#    plot(tma)
#    highlight(myFilter, recursive=TRUE)
# }
#
# \seealso{
#   See also the @see "RejectFilter" class.
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("GeneRejectFilter", function(input, layout, genes=NULL, ids=NULL, names=NULL, ...) {
  reject <- NULL;
  if (missing(input)) {
    input <- NULL;
    layout <- NULL;
  } else {
    if (missing(layout) || !inherits(layout, "Layout"))
      throw("Argument 'layout' must be of class Layout.");

    if (!is.null(ids)) {
      reject <- indexOf(layout, id=ids)
    } else if (!is.null(names)) {
      reject <- indexOf(layout, name=names);
    } else if (is.null(genes)) {
      throw("Either 'genes', 'ids' or 'names' must be given.");
    } else {
      geneGroups <- getGeneGroups(layout);
      if (is.null(geneGroups))
        throw("The specified layout has no genes specified.");
      reject <- getSpots(geneGroups, genes, unlist=TRUE);
    }
  }

  extend(RejectFilter(input=input, reject=reject, ...), "GeneRejectFilter", 
    layout = layout,
    genes  = genes
  )
})


setMethodS3("as.character", "GeneRejectFilter", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- data.class(this);
  s <- paste(s, ": Rejects a set of ", length(this$genes), " genes.", sep="");
  s;
})


setMethodS3("getIndex", "GeneRejectFilter", function(this) {
  incl <- getIndex(this$.input);
  idx <- setdiff(incl, this$reject);
  attr(idx, "max") <- attr(incl, "max");
  idx;
})



############################################################################
# HISTORY:
# 2002-02-26
# * Rewritten for setClassS3 and setMethodS3.
# 2001-07-19
# * Created!
############################################################################
