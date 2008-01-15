#########################################################################/**
# @RdocClass GeneAcceptFilter
#
# @title "The GeneAcceptFilter class"
#
# \description{
#  @classhierarchy
#
#   An GeneAcceptFilter asks its input for indices and let only those
#   indices through that corresponds to a given set of genes.
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
#   \item{...}{Any arguments accepted by the @see "AcceptFilter" constructor.}
#
#   Either 'genes', 'ids' or 'names' must be given.
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
#    # And among those only look at the duplicated genes.
#    genes <- getGeneGroups(layout)
#    duplicates <- which(getSizes(genes) == 2)
#    myFilter <- GeneAcceptFilter(fM, layout=layout, genes=duplicates, col="blue")
#    plot(tma)
#    highlight(myFilter, recursive=TRUE)
# }
#
# \seealso{
#   See also the @see "AcceptFilter" class.
# }
#
# \keyword{manip}
#*/#########################################################################
setConstructorS3("GeneAcceptFilter", function(input, layout, genes=NULL, ids=NULL, names=NULL, ...) {
  accept <- NULL;
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (missing(input)) {
    input <- NULL;
    layout <- NULL;
  } else {
    if (missing(layout) || !inherits(layout, "Layout"))
      throw("Argument 'layout' must be of class Layout.");
  
    if (!is.null(ids)) {
      genes <- indexOf(layout, id=ids)
    } else if (!is.null(names)) {
      genes <- indexOf(layout, name=names);
    } else if (is.null(genes)) {
      throw("Either 'genes', 'ids' or 'names' must be given.");
    } else {
      geneGroups <- getGeneGroups(layout);
      if (is.null(geneGroups))
        throw("The specified layout has no genes specified.");
      accept <- getSpots(geneGroups, genes, unlist=TRUE);
    }
  }

  extend(AcceptFilter(input=input, accept=accept, ...), "GeneAcceptFilter", 
    layout = layout,
    genes = genes
  )
})


setMethodS3("as.character", "GeneAcceptFilter", function(this) {
  s <- data.class(this);
  s <- paste(s, ": Accepts a set of ", length(this$genes), " genes.", sep="");
  s;
})



setMethodS3("getIndex", "GeneAcceptFilter", function(this) {
  incl <- getIndex(this$.input);
  idx <- intersect(incl, this$accept);
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
