# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
# MicroarrayArray
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
setConstructorS3("MicroarrayArray", function(array=NULL, layout=NULL) {
  if (!is.null(layout)) {
    if (!inherits(layout, "Layout"))
      throw("Argument 'layout' is not a Layout object: ", data.class(layout));
  }

  extend(Array(array), "MicroarrayArray", 
    .layout                   = layout,
    .saturated                = NULL,
    .view                     = c(AUTO.VIEW=0),
    AUTO.VIEW                 = c(AUTO.VIEW=0),
    DEFAULT.VIEW              = c(DEFAULT.VIEW=1),
    VECTOR.VIEW               = c(VECTOR.VIEW=2),
    SPOT.SLIDE.VIEW           = c(SPOT.SLIDE.VIEW=3),
    GENE.REPLICATE.SLIDE.VIEW = c(GENE.REPLICATE.SLIDE.VIEW=4),
    GENE.SLIDE.VIEW           = c(GENE.SLIDE.VIEW=5)
  )
})

setMethodS3("str", "MicroarrayArray", function(object, ...) {
  NextMethod();
})


setMethodS3("getLayout", "MicroarrayArray", function(this) {
  attr(this, ".layout");
})


setMethodS3("getView", "MicroarrayArray", function(this, call=NULL, set=FALSE) {
  view <- attr(this, ".view")[1];

  if (view == MicroarrayArray$AUTO.VIEW) {
    if (!is.null(call)) {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Automatically identify the view
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Check how many ,'s (commas) are used.
  #    call <- sys.call(sys.parent(n=1));
      # Exclude the function name and the object name
      call <- call[-c(1,2)];
      # Keep only the first three arguments that could be i, j and k.
      # Exclude the 'drop' arguments
      if (length(call) > 3)
  	call <- call[1:3];
      names <- names(call);
      if (!is.null(names))
  	call <- call[is.na(pmatch(names, "drop"))];
      len <- length(call);
      if (set == TRUE)
        len <- len - 1;
      # Note that len == 1 for both x[1] and x[], but the latter
      # has call[[1]] == ""...
      if (len == 1 && call[[1]] == "")
        len <- 0;
      if (len == 0) {
        view <- MicroarrayArray$DEFAULT.VIEW;
      } else if (len == 1) {
        view <- MicroarrayArray$VECTOR.VIEW;
      } else if (len == 2) {
  	view <- MicroarrayArray$SPOT.SLIDE.VIEW;
      } else if (len == 3) {
  	view <- MicroarrayArray$GENE.REPLICATE.SLIDE.VIEW;
      }
    } else {
      view <- MicroarrayArray$DEFAULT.VIEW;
    }
  }

  view;
})

setMethodS3("setView", "MicroarrayArray", function(this, newView) {
  if (is.null(newView))
    newView <- 0;
  if (!is.element(newView, 0:5))
    throw("Value of 'newView' is out of range: ", newView);
  attr(this, ".view") <- newView;
  invisible(this);
})

setMethodS3("pushView", "MicroarrayArray", function(this, view) {
  if (is.null(view))
    view <- 0;
  if (!is.element(view, 0:5))
    throw("Value of 'view' is out of range: ", view);
  attr(this, ".view") <- c(view, attr(this, ".view"));
  invisible(this);
})

setMethodS3("popView", "MicroarrayArray", function(this) {
  if (length(attr(this, ".view")) < 2)
    throw("Trying to pop from the view stack without pushing first.");
  view <- attr(this, ".view")[1];
  attr(this, ".view") <- attr(this, ".view")[-1];
  invisible(this);
})

setMethodS3("as.matrix", "MicroarrayArray", function(x) {
  # To please R CMD check...
  this <- x;

  unclass(getSpotOrGeneSlideValues(this))
})


setMethodS3("getSpotOrGeneSlideValues", "MicroarrayArray", abstract=TRUE);

setMethodS3("getSpotSlideValues", "MicroarrayArray", abstract=TRUE);
setMethodS3("setSpotSlideValues", "MicroarrayArray", abstract=TRUE);

setMethodS3("getGeneSlideValues", "MicroarrayArray", abstract=TRUE);
setMethodS3("setGeneSlideValues", "MicroarrayArray", abstract=TRUE);

setMethodS3("getGeneReplicateSlideValues", "MicroarrayArray", abstract=TRUE);
setMethodS3("setGeneReplicateSlideValues", "MicroarrayArray", abstract=TRUE);


setMethodS3("getSpotSlideValues", "ANY", function(x, spots=NULL, slides=NULL, drop=FALSE) {
  if (is.null(spots))
    spots <- seq(nrow(x))
  if (is.null(slides))
    slides <- seq(ncol(x))
  x[spots,slides, drop=drop];
})

setMethodS3("getGeneSlideValues", "ANY", function(x, genes=NULL, slides=NULL, drop=FALSE) {
  if (is.null(genes))
    genes <- seq(nrow(x))
  if (is.null(slides))
    slides <- seq(ncol(x))
  x[genes,slides, drop=drop];
})

setMethodS3("getGeneReplicateSlideValues", "ANY", function(x, genes=NULL, replicates=NULL, slides=NULL) {
  throw("Not implemented yet!")
})

setMethodS3("as.SpotSlideArray", "MicroarrayArray", function(this) {
  SpotSlideArray(getSpotSlideValues(this), layout=this$.layout);
})

setMethodS3("as.GeneSlideArray", "MicroarrayArray", function(this) {
  GeneSlideArray(getGeneSlideValues(this), layout=this$.layout);
})

setMethodS3("as.GeneReplicateSlideArray", "MicroarrayArray", function(this) {
  GeneReplicateSlideArray(getGeneReplicateSlideValues(this), layout=this$.layout);
})

setMethodS3("[", "MicroarrayArray", function(this, i=NULL, j=NULL, k=NULL, drop=TRUE) {
  view <- getView(this, call=sys.call(sys.parent(n=0)));
  if (view == MicroarrayArray$VECTOR.VIEW) {
    NextMethod();
  } else if (view == MicroarrayArray$SPOT.SLIDE.VIEW) {
    getSpotSlideValues(this, spots=i, slides=j, drop=drop);
  } else if (view == MicroarrayArray$GENE.REPLICATE.SLIDE.VIEW) {
    getGeneReplicateSlideValues(this, genes=i, replicates=j, slides=k, drop=drop);
  } else if (view == MicroarrayArray$GENE.SLIDE.VIEW) {
    getGeneSlideValues(this, genes=i, slides=j, drop=drop);
  }
}, createGeneric=FALSE)


setMethodS3("[<-", "MicroarrayArray", function(this, i=NULL, j=NULL, k=NULL, value) {
  view <- getView(this, call=sys.call(sys.parent(n=0)), set=TRUE);
  if (view == MicroarrayArray$VECTOR.VIEW) {
    NextMethod();
  } else if (view == MicroarrayArray$SPOT.SLIDE.VIEW) {
    setSpotSlideValues(this, spots=i, slides=j, value=value);
  } else if (view == MicroarrayArray$GENE.REPLICATE.SLIDE.VIEW) {
    setGeneReplicateSlideValues(this, genes=i, replicates=j, slides=k, value=value);
  } else if (view == MicroarrayArray$GENE.SLIDE.VIEW) {
    setGeneSlideValues(this, genes=i, slides=j, value=value);
  }
}, createGeneric=FALSE);


# Include these for backward compatibility.
setMethodS3("pushView", "ANY", function(this, view) {
  invisible(this);
})

setMethodS3("popView", "ANY", function(this) {
  invisible(this);
})


setMethodS3("getView", "ANY", function(this) {
  MicroarrayArray$SPOT.SLIDE.VIEW;
})



############################################################################
# HISTORY:
# 2002-12-31
# o Added drop=FALSE arguments where applicable.
# 2002-11-12
# o Added support for retrieving and assigning as vectors, e.g. 
#   x[1:10] <- 0. Updated with a new view VECTOR.VIEW.
# o All view values are now named to simplifying debugging etc.
# o Update the getView() function with the argument 'set=FALSE'. Any 'set' 
#   functions should call getView() with 'set=TRUE', because they have one
#   more argument, namely the value.
# o Added abstract methods of setXXXSlideArray().
# 2002-11-07
# o Added pushView() and popView().
# o Added abstract "["() and "[<-"().
# 2002-11-03
# o Created!
############################################################################
