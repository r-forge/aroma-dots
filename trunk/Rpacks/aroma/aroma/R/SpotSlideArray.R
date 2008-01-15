# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# SpotSlideArray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
setConstructorS3("SpotSlideArray", function(array=NULL, layout=NULL) {
  if (!is.null(array))
    array <- as.matrix(array);

  extend(MicroarrayArray(array=array, layout=layout), "SpotSlideArray"
  )
})



setMethodS3("getView", "SpotSlideArray", function(this, ...) {
  view <- NextMethod("getView", this, ...);
  if (view == this$DEFAULT.VIEW)
    view <- this$SPOT.SLIDE.VIEW;
  view;
})


setMethodS3("as.matrix", "SpotSlideArray", function(x) {
  # To please R CMD check...
  this <- x;

  unclass(this)
})


setMethodS3("getSpotOrGeneSlideValues", "SpotSlideArray", function(this, ...) {
  getSpotSlideValues(this, ...);
})



setMethodS3("getSpotSlideValues", "SpotSlideArray", function(this, spots=NULL, slides=NULL, drop=FALSE) {
  uthis <- unclass(this);

  if (is.null(spots))
    spots <- seq(nrow(uthis));

  if (is.null(slides))
    slides <- seq(ncol(uthis));

  uthis[spots, slides, drop=drop];
})


setMethodS3("setSpotSlideValues", "SpotSlideArray", function(this, spots=NULL, slides=NULL, value) {
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




setMethodS3("getGeneReplicateSlideValues", "SpotSlideArray", function(this, genes=NULL, replicates=NULL, slides=NULL, drop=FALSE) {
  uthis <- unclass(this);   # Calling unclass() takes time!

  layout <- getLayout(this);
  geneGroups <- getGeneGroups(layout);

  if (is.null(genes))
    genes <- seq(geneGroups);

  # Get a list of element each containing the indices for each gene.
  geneSpots <- getSpots(geneGroups, groups=genes);

  # Keep only those replicates that are specified by j. Note that it
  # is possible to specify too many replicates.
  if (!is.null(replicates)) {
    sizes <- unlist(lapply(geneSpots, FUN=length));
    for (l in seq(along=geneSpots)) {
      include <- intersect(seq(sizes[l]), replicates);
      spots[[l]] <- spots[[l]][include];
    }
  }

  # Get the slides to be included
  if (is.null(slides))
    slides <- seq(ncol(uthis));

  nbrOfSlides <- length(slides);

  array <- list();
  for (l in seq(along=geneSpots))
    array[[l]] <- matrix(uthis[geneSpots[[l]],slides], ncol=nbrOfSlides);

  array;
})


setMethodS3("getGeneSlideValues", "SpotSlideArray", function(this, genes=NULL, slides=NULL, drop=FALSE) {
  throw("Not implement yet!");
})




############################################################################
# HISTORY:
# 2002-12-31
# o Added drop=FALSE arguments where applicable.
# 2002-11-12
# o Update the getView() function to pass on *any* extra arguments to 
#   getView() in class MicroarrayArray.
# o Added setSpotSlideArray().
# o Added getSpotOrGeneSlideArray().
# 2002-11-06
# o Added as.GeneSlideArray().
# 2002-11-01
# o Created!
############################################################################

