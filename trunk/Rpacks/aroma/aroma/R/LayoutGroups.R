#########################################################################/**
# @RdocClass LayoutGroups
#
# @title "The LayoutGroups class"
#
# \description{
#  @classhierarchy
#
#  The LayoutGroups class is a companion class for the Layout class.
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @examples "../incl/LayoutGroups.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setConstructorS3("LayoutGroups", function(layout=NULL) {
  if (!is.null(layout) & !inherits(layout, "Layout"))
    throw("Argument 'layout' must of class Layout: ", data.class(layout));
    
  extend(Object(), "LayoutGroups",
    layout=layout
  )
}, abstract=TRUE);


setMethodS3("as.character", "LayoutGroups", function(object) {
  s <- paste(data.class(object), ": ", nbrOfGroups(object),
             " groups defined.", sep="");
  s;
})

setMethodS3("getLayout", "LayoutGroups", function(object) {
  object$layout;
})

setMethodS3("setLayout", "LayoutGroups", function(object, layout) {
  if (!inherits(layout, "Layout"))
    throw("Argument 'layout' must of class Layout: ", data.class(layout));
  object$layout <- layout;
})

setMethodS3("getNames", "LayoutGroups", function(object) {
  NULL;
}, abstract=TRUE)

setMethodS3("nbrOfGroups", "LayoutGroups", function(object) {
  length(getSpots(object));
})

setMethodS3("seq", "LayoutGroups", function(object) {
  1:nbrOfGroups(object);
})

setMethodS3("getSpots", "LayoutGroups", function(object, groups=NULL, unlist=FALSE) {
  NULL;
}, abstract=TRUE)


setMethodS3("getGroupValues", "LayoutGroups", function(object, data, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) groups <- seq(object);
  data <- as.matrix(data);
  spots <- getSpots(object, groups=groups);
  l <- list();
  for (k in 1:length(spots)) {
    l[[k]] <- data[spots[[k]],];
  }
  names(l) <- getNames(object)[groups];
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})

setMethodS3("getSpotValues", "LayoutGroups", function(object, data, groups=NULL, unlist=FALSE) {
  res <- rep(NA, nbrOfSpots(getLayout(object)));
  spots <- getSpots(object, groups=groups);
  for (k in 1:length(spots)) {
    res[spots[[k]]] <- data[k];
  }
  res;
})


setMethodS3("apply", "LayoutGroups", function(object, data, FUN, groups=NULL, unlist=FALSE, ...) {
  if (is.null(groups)) groups <- seq(object);
  values <- getGroupValues(object, data, groups=groups);
  l <- lapply(values, FUN=FUN, ...);
  names(l) <- getNames(object)[groups];
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
}, abstract=TRUE)





#########################################################################/**
# @RdocClass PlateGroups
#
# @title "The PlateGroups class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
#  \item{plateDef}{
#   If a scalar (@numeric of length one), each plate group is assumed to
#   contain that number of spots/clones.
#   Can also be a @data.frame.
#   If @NULL, the plates are obtained from the Layout 
#   object, which then is required to contain plate information.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @examples "../incl/LayoutGroups.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setConstructorS3("PlateGroups", function(layout=NULL, plateDef=NULL) {
  groups <- NULL;
  if (!is.null(layout)) {
    if (is.null(plateDef)) {
      if (!hasPlates(layout))
        throw("If argument 'plateDef' is not defined, then argument 'layout' must have plates defined, which is not the case.");
      # Assume all printtips have the same plate layout as printtip #1
      plateNbr <- getPlateNumber(layout)[1:(gridSize(layout)+1)];
      idx <- which(diff(plateNbr) != 0);
      last <- idx * nbrOfGrids(layout);
      names(last) <- getPlate(layout, idx);
      plateDef <- last;
    }
    
    if (is.numeric(plateDef)) {
      if (length(plateDef) == 1) {
        plateSize <- plateDef
        nbrOfPlates <- nbrOfSpots(layout) / plateSize
        first <- (0:(nbrOfPlates-1))*plateSize + 1
        last <- first + plateSize - 1
        names <- 1:length(first)
      } else {
        names <- names(plateDef)
        if (plateDef[1] == 1) {
          first <- unique(plateDef)
          last <- c(first[-1]-1, nbrOfSpots(layout))
        } else {
          last <- unique(plateDef)
          first <- c(1, last[-length(last)]+1)
        }
      }
      mat <- t(matrix(c(first, last), ncol=2))
      groups <- as.list(as.data.frame(mat, optional=TRUE))
      names(groups) <- names;
    } else if (is.data.frame(plateDef)) {
      groups <- as.list(plateDef);
    }
  }
  
  extend(LayoutGroups(layout), "PlateGroups",
    groups=groups
  )
})

setMethodS3("as.character", "PlateGroups", function(object) {
  s <- NextMethod("as.character");
  first <- getFirst(object);
  last <- getLast(object);
  names <- getNames(object);
  if (!is.null(names))
    names <- paste("'", names, "' ", sep="");
  tmp <- paste(names ,"(", first, "-", last, ")", sep="", collapse=", ");
  s <- paste(s, " The plates groups are ", tmp, ".", sep="");
  s;
})

setMethodS3("getFirst", "PlateGroups", function(object) {
  matrix(unlist(object$groups), nrow=2)[1,]
})

setMethodS3("getLast", "PlateGroups", function(object) {
  matrix(unlist(object$groups), nrow=2)[2,]
})

setMethodS3("getNames", "PlateGroups", function(object) {
  names(object$groups);
})

setMethodS3("nbrOfGroups", "PlateGroups", function(object) {
  length(object$groups)
})

setMethodS3("getSizes", "PlateGroups", function(object) {
  unlist(lapply(object$groups, FUN=length))
})


setMethodS3("getPrintorderIndices", "PlateGroups", function(object, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(object))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(object)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(object))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the plate names.");
    groups <- match[ok];
  }

  plates <- object$groups[groups]
  l <- lapply(plates, FUN=function(x) x[1]:x[2])
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})


setMethodS3("getSpots", "PlateGroups", function(object, groups=NULL, unlist=FALSE) {
  l0 <- getPrintorderIndices(object, groups=groups, unlist=unlist)
  m <- as.vector(toPrintorderMatrix(object$layout))
  l <- lapply(l0, FUN=function(i) m[i])
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})







############################################################################
# SuperGroups
############################################################################
setConstructorS3("SuperGroups", function(layoutGroups=NULL, groups=NULL) {
  layout <- NULL;
  if (!is.null(layoutGroups)) {
    if (!inherits(layoutGroups, "LayoutGroups"))
      throw("Argument 'layoutGroups' is not of class LayoutGroups.")
    layout <- getLayout(layoutGroups);
    if (is.null(groups))
      groups <- as.list(1:nbrOfGroups(layoutGroups));
  }

  groups <- as.list(groups);
  
  extend(LayoutGroups(layout), "SuperGroups",
    layoutGroups=layoutGroups,
    groups=groups
  )
});


setMethodS3("as.character", "SuperGroups", function(object) {
  s <- paste(data.class(object), ": ", nbrOfGroups(object),
             " groups defined.", sep="");
  s;
})

setMethodS3("getNames", "SuperGroups", function(object) {
  NULL;
}, abstract=TRUE)

setMethodS3("getGroups", "SuperGroups", function(object) {
  object$groups;
}, abstract=TRUE)

setMethodS3("setGroups", "SuperGroups", function(object, groups) {
  object$groups <- as.list(groups);
}, abstract=TRUE)

setMethodS3("nbrOfGroups", "SuperGroups", function(object) {
  length(object$groups);
}, abstract=TRUE)


setMethodS3("getSizes", "SuperGroups", function(object) {
  unlist(lapply(getSpots(object), FUN=length))
})


setMethodS3("getSpots", "SuperGroups", function(object, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(object))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(object)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else {
    throw("Argument 'groups' must be numeric.");
  }

  layoutGroups <- object$layoutGroups;
  groups <- object$groups[groups];
  l <- list();
  for (k in 1:length(groups)) {
    cluster <- groups[[k]];
    spots <- getSpots(layoutGroups, cluster, unlist=TRUE)
    l[[k]] <- spots;
  }
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
}, abstract=TRUE)




#########################################################################/**
# @RdocClass SlideRowGroups
#
# @title "The SlideRowGroups class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
#  \item{groups}{A @list of length equal to the number of rows of spots
#   on each slide and each element contains the spot indices in each row.
#   If @NULL, the rows are obtained from the Layout object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @examples "../incl/LayoutGroups.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setConstructorS3("SlideRowGroups", function(layout=NULL, groups=NULL) {
  if (!is.null(layout) && is.null(groups)) {
    xy <- toXYMatrix(layout);
    groups <- list();
    for (k in 1:nrow(xy))
      groups <- c(groups, list(xy[k,]));
  }

  extend(LayoutGroups(layout), "SlideRowGroups",
    groups=groups
  )
})

setMethodS3("as.character", "SlideRowGroups", function(object) {
  s <- NextMethod("as.character");
  idx <- getSpots(object);
  first <- unlist(lapply(idx, FUN=function(x) x[1]));
  last <- unlist(lapply(idx, FUN=function(x) x[length(x)]));
  tmp <- paste(first, last, sep="-", collapse=", ");
  s <- paste(s, " The row groups are ", tmp, ".", sep="");
  s;
})

setMethodS3("getNames", "SlideRowGroups", function(object) {
  names(object$groups);
})

setMethodS3("nbrOfGroups", "SlideRowGroups", function(object) {
  length(object$groups)
})

setMethodS3("getSizes", "SlideRowGroups", function(object) {
  unlist(lapply(object$groups, FUN=length))
})


setMethodS3("getSpots", "SlideRowGroups", function(object, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(object))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(object)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(object))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the row names.");
    groups <- match[ok];
  }
  
  l <- object$groups;
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})






#########################################################################/**
# @RdocClass SlideColumnGroups
#
# @title "The SlideColumnGroups class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#  \item{layout}{A @see "Layout" object specifying the layout of a set of
#   microarrays.}
#  \item{groups}{A @list of length equal to the number of rows of spots
#   on each slide and each element contains the spot indices in each row.
#   If @NULL, the rows are obtained from the Layout object.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
# @author
#
# @examples "../incl/LayoutGroups.Rex"
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setConstructorS3("SlideColumnGroups", function(layout=NULL, groups=NULL) {
  if (!is.null(layout) && is.null(groups)) {
    xy <- toXYMatrix(layout);
    groups <- list();
    for (k in 1:ncol(xy))
      groups <- c(groups, list(xy[,k]));
  }

  extend(LayoutGroups(layout), "SlideColumnGroups",
    groups=groups
  )
})

setMethodS3("as.character", "SlideColumnGroups", function(object) {
  s <- NextMethod("as.character");
  idx <- getSpots(object);
  first <- unlist(lapply(idx, FUN=function(x) x[1]));
  last <- unlist(lapply(idx, FUN=function(x) x[length(x)]));
  tmp <- paste(first, last, sep="-", collapse=", ");
  s <- paste(s, " The column groups are ", tmp, ".", sep="");
  s;
})

setMethodS3("getNames", "SlideColumnGroups", function(object) {
  names(object$groups);
})

setMethodS3("nbrOfGroups", "SlideColumnGroups", function(object) {
  length(object$groups)
})

setMethodS3("getSizes", "SlideColumnGroups", function(object) {
  unlist(lapply(object$groups, FUN=length))
})


setMethodS3("getSpots", "SlideColumnGroups", function(object, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(object))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(object)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(object))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the column names.");
    groups <- match[ok];
  }
  
  l <- object$groups;
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})







############################################################################
# HISTORY:
# 2003-09-19
# o Extracted PrintdipGroups and PrinttipGroups into its own source file.
# 2003-05-04
# o Updated the Rdoc's with argument specifications etc.
# 2003-04-08
# o Removed all missing links in the Rdoc comments.
# 2002-12-05
# o Added SlideRowGroups and SlideColumnGroups.
# 2002-11-27
# o Added the possibility to create a GeneGroups object where the within-
#   slide replication specification overrides the same information given
#   by the Layout object. Added the (obvious) specification "none".
# 2002-06-24
# o Example fix: Updated several of the example codes.
# 2002-05-10
# o Added seq(), getGroupValues() and getSpotValues() to LayoutGroups.
# o Added SuperGroups() to create clusters of LayoutGroups.
# 2002-05-04
# o Added PrintdipGroups()!
# 2002-05-03
# o Added support for specifying genes by there names.
# o Renamed ReplicateGroups() to GeneGroups().
# o getSizes() now returnes a vector as expected (not a list).
# o Renames getSpotIndices() to getSpots().
# o If a ReplicateGroups object is created an no replicate specification is
#   found, then no replicates are assumed and a warning is thrown.
# o Added getPrintorderIndices() to PlateGroups.
# 2002-05-02
# o Added getNames() and getSizes().
# o Added some Rdoc comments, mostly to get something and a running example.
# o Added PrinttipGroups, PlateGroups and ReplicateGroups. What other
#   LayoutGroups are there?
# o Created. With the introduction of plate groups also it feels like one
#   has to generalize the concept of groups; printtip groups, plate groups,
#   replicate groups. Still don't know how to introduce the concept of
#   slides; there is no information about the number of slides in Layout.
############################################################################
