# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# GeneReplicateSlideArray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# In a GeneReplicateSlideArray the microarray data is stored as a list
# where each element contains a matrix where each row represents a
# replicate and each column a slide.
setConstructorS3("GeneReplicateSlideArray", function(array=NULL, layout=NULL) {
  if (!is.null(array)) {
    if (!is.list(array)) {
      array <- as.list(array);
      array <- lapply(array, FUN=as.matrix);
    }
  }

  extend(MicroarrayArray(array=array, layout=layout), "GeneReplicateSlideArray"
  )
})

setMethodS3("ncol", "GeneReplicateSlideArray", function(this) {
  ncol(unclass(this)[[1]]);
})

setMethodS3("nrow", "GeneReplicateSlideArray", function(this) {
  nbrOfSpots(this$.layout);
})



setMethodS3("getView", "GeneReplicateSlideArray", function(this, ...) {
  view <- NextMethod("getView", this, ...);
  if (view == this$DEFAULT.VIEW)
    view <- this$GENE.REPLICATE.SLIDE.VIEW;
  view;
})


setMethodS3("getSpotOrGeneSlideValues", "GeneReplicateSlideArray", function(this, ...) {
  getGeneSlideValues(this, ...);
})


setMethodS3("getSpotSlideValues", "GeneReplicateSlideArray", function(this, spots=NULL, slides=NULL, drop=FALSE) {
  layout <- getLayout(this);
  uthis <- unclass(this);

  if (is.null(spots))
    spots <- seq(nbrOfSpots(layout));

  if (is.null(slides))
    slides <- seq(ncol(uthis[[1]]));

  # For each gene, get the spot indices for that gene.
  genes <- getGeneGroups(layout);
  spots <- getSpots(genes);
  
  map <- rep(NA, nbrOfSpots(layout));
  for (k in 1:length(spots))
    map[spots[[k]]] <- k;

  nbrOfSpots  <- length(spots);
  nbrOfSlides <- length(slides);

  array <- matrix(NA, nrow=nbrOfSpots, ncol=nbrOfSlides);
  rownames(array) <- spots;
  colnames(array) <- slides;

  for (l in seq(along=spots)) {
    spot      <- spots[l];
    gene      <- map[spot];
    replicate <- sum(map[1:spot] == gene);
    value     <- uthis[[gene]][,slides, drop=FALSE];
    replicate.max <- max(replicate);
    if (nrow(value) < replicate.max) {
      value     <- rep(value, length.out=replicate.max*nbrOfSlides);
      value     <- matrix(value, ncol=nbrOfSlides);
    }
    array[l,] <- value[replicate,];
  }
  array;
})



setMethodS3("getGeneReplicateSlideValues", "GeneReplicateSlideArray", function(this, genes=NULL, replicates=NULL, slides=NULL, drop=FALSE) {
  if (is.null(genes))
    genes <- seq(length(this));

  uthis <- unclass(this);
  array <- uthis[genes];


  # Keep only those replicates that are specified by j. Note that it
  # is possible to specify too many replicates.
  if (!is.null(replicates)) {
    sizes <- unlist(lapply(array, FUN=nrow));
    for (l in seq(along=array)) {
      include <- intersect(seq(sizes[l]), replicates);
      array[[l]] <- array[[l]][include,,drop=FALSE];
    }
  }

  # Get the slides to be included
  if (!is.null(slides)) {
    for (l in seq(along=array))
      array[[l]] <- array[[l]][,slides,drop=FALSE];
  }
  array;
})


setMethodS3("getGeneSlideValues", "GeneReplicateSlideArray", function(this, genes=NULL, slides=NULL, drop=FALSE) {
  uthis <- unclass(this);

  if (is.null(genes))
    genes <- seq(nrow(uthis));

  if (is.null(slides))
    slides <- seq(ncol(uthis));

  array <- matrix(NA, nrow=length(genes), ncol=length(slides));
  rownames(array) <- genes;
  colnames(array) <- slides;
  for (k in seq(along=genes)) {
    array[k,] <- uthis[[k]][1,slides];
  }
  array;
})



############################################################################
# HISTORY:
# 2002-11-12
# o Update the getView() function to pass on *any* extra arguments to 
#   getView() in class MicroarrayArray.
# o Added getSpotOrGeneSlideArray().
# 2002-11-06
# o Added as.GeneSlideArray().
# 2002-11-01
# o Created!
############################################################################
