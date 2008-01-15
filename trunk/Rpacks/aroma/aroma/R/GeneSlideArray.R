# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GeneSlideArray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# In a GeneSlideArray the microarray data is stored as a list
# where each element contains a matrix where each row represents a
# replicate and each column a slide.
setConstructorS3("GeneSlideArray", function(array=NULL, layout=NULL) {
  if (!is.null(array)) {
    if (!is.matrix(array))
      array <- as.matrix(array);
  }

  extend(MicroarrayArray(array=array, layout=layout), "GeneSlideArray"
  )
})

setMethodS3("ncol", "GeneSlideArray", function(this) {
  ncol(unclass(this));
})

setMethodS3("nrow", "GeneSlideArray", function(this) {
  nbrOfSpots(this$.layout);
})


setMethodS3("getView", "GeneSlideArray", function(this, ...) {
  view <- NextMethod("getView", this, ...);
  if (view == this$DEFAULT.VIEW)
    view <- this$GENE.SLIDE.VIEW;
  view;
})

setMethodS3("as.matrix", "GeneSlideArray", function(x) {
  # To please R CMD check...
  this <- x;

  unclass(this)
})


setMethodS3("getSpotOrGeneSlideValues", "GeneSlideArray", function(this, ...) {
  getGeneSlideValues(this, ...);
})

setMethodS3("getSpotSlideValues", "GeneSlideArray", function(this, spots=NULL, slides=NULL) {
  layout <- getLayout(this);
  uthis <- unclass(this);

  if (is.null(spots))
    spots <- seq(nbrOfSpots(layout));

  if (is.null(slides))
    slides <- seq(ncol(uthis))
  else
    uthis <- uthis[,slides, drop=FALSE];
  
  # For each gene, get the spot indices for that gene.
  geneGroups <- getGeneGroups(layout);
  geneSpots  <- getSpots(geneGroups);

  res <- matrix(NA, nrow=nbrOfSpots(layout), ncol=length(slides));
  for (k in seq(along=geneSpots)) {
    value <- uthis[k,];
    for (spot in geneSpots[[k]])
      res[spot,] <- value;
  }

  res[spots, , drop=FALSE];
})



setMethodS3("getGeneReplicateSlideValues", "GeneSlideArray", function(this, genes=NULL, replicates=NULL, slides=NULL) {
  if (is.null(genes))
    genes <- seq(length(this));

  uthis <- unclass(this);
  array <- uthis[genes, , drop=FALSE];

  if (is.null(replicates))
    replicates <- 1;

  if (is.null(slides))
    slides <- seq(ncol(array));

  nbrOfReplicates <- length(replicates);
  nbrOfSlides     <- length(slides);

  res <- list();
  for (l in seq(along=genes)) {
    res[[l]] <- matrix(array[l,slides], 
                       nrow=nbrOfReplicates, ncol=nbrOfSlides, byrow=TRUE);
  }

  res;
})


setMethodS3("getGeneSlideValues", "GeneSlideArray", function(this, genes=NULL, slides=NULL) {
  uthis <- unclass(this);

  if (is.null(genes))
    genes <- seq(nrow(uthis));

  if (is.null(slides))
    slides <- seq(ncol(uthis));

  uthis[genes, slides, drop=FALSE];
})


setMethodS3("setGeneSlideValues", "GeneSlideArray", function(this, spots=NULL, slides=NULL, value) {
  uthis <- unclass(this);

  if (is.null(spots))
    spots <- seq(nrow(uthis));

  if (is.null(slides))
    slides <- seq(ncol(uthis));

  uthis[spots, slides] <- value;
  for (name in names(attributes(this)))
    attr(uthis, name) <- attr(this, name);
  uthis;
})




############################################################################
# HISTORY:
# 2002-11-12
# o Update the getView() function to pass on *any* extra arguments to 
#   getView() in class MicroarrayArray.
# o Added setGeneSlideArray().
# o Added getSpotOrGeneSlideArray().
# 2002-11-06
# o Added support for view GENE.SLIDE.VIEW.
# 2002-11-05
# o Created!
############################################################################
