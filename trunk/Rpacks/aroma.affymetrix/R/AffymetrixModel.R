###########################################################################/**
# @RdocClass AffymetrixModel
#
# @title "The AffymetrixModel class"
#
# \description{
#  @classhierarchy
#
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataset}{}
#   \item{path}{}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# \seealso{
# }
#
# @author
#*/###########################################################################
setConstructorS3("AffymetrixModel", function(dataset=NULL, path=NULL, ...) {
  if (!is.null(dataset)) {
    if (!inherits(dataset, "AffymetrixDataset"))
      throw("Argument 'dataset' is not an AffymetrixDataset object: ", 
                                                           class(dataset)[1]);
  }

  extend(Object(), "AffymetrixModel", 
    .dataset = dataset,
    .path = path,
    .fit = NULL
  )
}, abstract=TRUE)


setMethodS3("as.character", "AffymetrixModel", function(this, ...) {
  s <- paste(class(this)[1], ":", sep="");
  s <- paste(s, " ", as.character(this$.dataset), sep="");
  s <- paste(s, " ", as.character(this$.fit), sep="");
  s;
})

setMethodS3("getDataset", "AffymetrixModel", function(this, ...) {
  this$.dataset;
})

setMethodS3("nbrOfArrays", "AffymetrixModel", function(this, ...) {
  nbrOfArrays(this$.dataset);
})

setMethodS3("nbrOfProbes", "AffymetrixModel", function(this, ...) {
  nbrOfProbes(this$.dataset);
})

setMethodS3("nbrOfProbesets", "AffymetrixModel", function(this, ...) {
  nbrOfProbesets(this$.dataset);
})

setMethodS3("getPath", "AffymetrixModel", function(this, ...) {
  this$.path;
})

setMethodS3("getFit", "AffymetrixModel", abstract=TRUE);

setMethodS3("fit", "AffymetrixModel", abstract=TRUE);

setMethodS3("fromAffymetrixDataset", "AffymetrixModel", abstract=TRUE);

setMethodS3("findUnitsDone", "AffymetrixModel", abstract=TRUE);


############################################################################
# HISTORY:
# 2006-04-05
# o Created from RlmmModel.R.
############################################################################

