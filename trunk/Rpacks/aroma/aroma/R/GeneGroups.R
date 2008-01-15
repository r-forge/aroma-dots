#########################################################################/**
# @RdocClass GeneGroups
#
# @title "The GeneGroups class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{layout}{A Layout object.}
#   \item{specification}{(optional) @character string specifying what type of 
#     within-slide replication the slides are printed with. If @NULL,
#     the within-slide replication will be extracted from the names or the
#     ids of the Layout object. Possible values are 
#     \code{"none"} (no within-slide replications), 
#     \code{"neighboring-pairs"} (duplicated horizontal neighboring pairs) and
#     \code{"top-bottom"} (duplicated pairs over replicated tiles).
#   }
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
#*/#########################################################################
setConstructorS3("GeneGroups", function(layout=NULL, specification=NULL) {
  groups <- NULL;
  if (!is.null(layout)) {
    nbrOfSpots <- nbrOfSpots(layout);

    if (!is.null(specification)) {
      # If specification is specified, use that information...
      if (is.character(specification)) {
  	type <- specification;
  	type <- tolower(type);
  	if (type == "none") {
          groups <- as.list(1:nbrOfSpots);
          names(groups) <- 1:nbrOfSpots;
        } else if (type == "neighboring-pairs" || type == "duplicates") {
  	  mat <- matrix(1:nbrOfSpots, nrow=2, byrow=FALSE);
  	  df <- as.data.frame(mat, optional=TRUE);
  	  groups <- as.list(df);
  	  names(groups) <- 1:(nbrOfSpots %/% 2);
  	} else if (type == "top-bottom") {
  	  mat <- matrix(1:nbrOfSpots, ncol=2, byrow=FALSE);
  	  df <- as.data.frame(mat, optional=TRUE);
  	  groups <- as.list(df);
  	  names(groups) <- 1:(nbrOfSpots %/% 2);
  	} else {
  	  throw("Argument 'specification' contains an unknown type string: ", type);
  	}
      } else if (is.list(specification)) {
        group <- specification;
      } else {
        throw("Unknown value of argument 'specification'. Must be a type string or a list.");
      }
    } else {
      # Try to extract the replicate information from the Layout object
      if (!hasNames(layout) && !hasIDs(layout)) {
        warning("If argument 'specification' is not defined, then argument 'layout' must have names or ids defined, other no within-slide replicates are assumed.");
        groups <- as.list(1:nbrOfSpots);
        names(groups) <- 1:nbrOfSpots;
      } else {
        x1 <- x2 <- u1 <- u2 <- NULL;
        if (hasNames(layout)) {
          x1 <- getName(layout);
          u1 <- unique(x1);
        }
        
        if (hasIds(layout)) {
          x2 <- getId(layout);
          u2 <- unique(x2);
        }
  
        n1 <- length(u1)
        n2 <- length(u2)
        if (n1 > n2) {
          u <- u1
          x <- x1
          n <- n1
        } else {
          u <- u2
          x <- x2
          n <- n2
        }
        rm(u1,u2,x1,x2,n1,n2)
        m <- match(x, u);
        groups <- list(); groups[[n+1]] <- NA;
        ks <- which(!is.na(m))
        for(k in ks) {
          i <- m[k];
          groups[[i]] <- c(groups[[i]], k);
        }
        groups[[n+1]] <- NULL;
        names(groups) <- u;
      }
    } # if (is.null(specification))
  } # if (!is.null(layout))
  
  extend(LayoutGroups(layout), "GeneGroups",
    groups=groups
  )
})


setMethodS3("as.character", "GeneGroups", function(this) {
  s <- NextMethod("as.character");
  idx <- getSpots(this);
  len <- unlist(lapply(idx, FUN=length));
  len <- sort(len);
  first <- c(1, which(c(diff(len)) != 0)+1, length(idx)+1);
  len <- len[first[-length(first)]];
  count <- diff(first);
  explicit <- which(count == 1);
  names <- rep("", length(len));
  names[explicit] <- paste("(", names(len)[explicit], ") ", sep="");
  tmp <- paste(count, " genes ", names, "with ", len, " replicates", sep="", collapse=", ");
  s <- paste(s, " ", tmp, ".", sep="");
  s;
})



#########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the unique gene names"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns a @vector of length \code{nbrOfGroups(this)}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getNames", "GeneGroups", function(this) {
  names(this$groups);
})



#########################################################################/**
# @RdocMethod nbrOfGroups
#
# @title "Gets the number of unique genes"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("nbrOfGroups", "GeneGroups", function(this) {
  length(this$groups)
})


#########################################################################/**
# @RdocMethod getSizes
#
# @title "Gets the number of replicates for each gene"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns an @integer @vector of length \code{nbrOfGroups(this)}.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSizes", "GeneGroups", function(this) {
  unlist(lapply(this$groups, FUN=length))
})




#########################################################################/**
# @RdocMethod getSpots
#
# @title "Gets a list of items each containing spot indices for a group"
#
# @synopsis
#
# \description{
#  Gets a list of items each containing spot indices for a group.
# }
#
# \arguments{
#   \item{groups}{An optional @vector of group indices specifying for
#     which groups the spot indicies should be returned.
#     If @NULL, all groups are used.}
# }
#
# \value{
#   Returns a @list of length \code{nbrOfGroups(this)}.
# }
#
# @author
#
# @examples "../incl/GeneGroups.getSpots.Rex"
#
# \seealso{
#   @seemethod "indexOf", @seemethod "setId".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getSpots", "GeneGroups", function(this, groups=NULL, unlist=FALSE) {
  if (is.null(groups)) {
    groups <- seq(nbrOfGroups(this))
  } else if (is.numeric(groups)) {
    if (any(groups < 1 || groups > nbrOfGroups(this)))
      throw("Argument 'groups' contain a value that is out of range.")
  } else if (is.character(groups)) {
    match <- match(unique(groups), getNames(this))
    ok <- !is.na(match)
    if (sum(!ok) != 0)
      warning("Some of the names asked for where not found among the gene names.");
    groups <- match[ok];
  }

  l <- this$groups[groups]
  if (unlist == TRUE) unlist(l,use.names=FALSE) else l
})




#########################################################################/**
# @RdocMethod getReplicates
#
# @title "Gets all spots that exist with a certain number of replicates"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \arguments{
#   \item{nbrOfReplicates}{The number of replicates to be selected for, i.e.
#     the returned genes all have this number of replicates.}
#   \item{value}{A @vector containing the values to be returned ordered in
#     the same way as the spot indicies. If @NULL, the spot indices are
#     returned.}
#   \item{asMatrix}{If @TRUE, the returned values are returned as a
#     matrix with \code{nbrOfReplicates} columns and where each row 
#     represents one gene and the rows are named as the genes, 
#     otherwise a @list is returned where each element represents one gene
#     and contains \code{nbrOfReplicates} values and the elements are named 
#     as the genes.}
# }
#
# \value{
#   Returns a @list or a @matrix.
# }
#
# @author
#
# \seealso{
#   Shorthand versions: @seemethod "getDuplicates", 
#   @seemethod "getTriplicates", @seemethod "getQuadruplicates".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getReplicates", "GeneGroups", function(this, nbrOfReplicates, value=NULL, asMatrix=TRUE) {
  idx <- getSpots(this);

  # Get all spots of length 'nbrOfReplicates'
  len <- unlist(lapply(idx, FUN=length));
  tups <- idx[(len == nbrOfReplicates)];

  # Return spot indices or values.
  if (!is.null(value))
    tups <- lapply(tups, FUN=function(idx) value[idx]);

  # Return as a list or as a matrix.
  if (asMatrix) {
    rownames <- names(tups);
    tups <- matrix(unlist(tups), ncol=nbrOfReplicates, byrow=TRUE);
    rownames(tups) <- rownames;
  }

  tups;
}) # getReplicates()




#########################################################################/**
# @RdocMethod getDuplicates
#
# @title "Gets all duplicated spots"
#
# @synopsis
#
# \description{
#  @get "title". This is a shortcut for calling
#  @seemethod "getReplicates" with \code{nbrOfReplicates=2}.
# }
#
# \arguments{
#   \item{...}{Arguments accepted by @seemethod "getReplicates".}
# }
#
# \value{
#   Returns a @list or a @matrix.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReplicates".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getDuplicates", "GeneGroups", function(this, ...) {
  getReplicates(this, nbrOfReplicates=2, ...);
}) # getDuplicates()



#########################################################################/**
# @RdocMethod getTriplicates
#
# @title "Gets all triplicated spots"
#
# @synopsis
#
# \description{
#  @get "title". This is a shortcut for calling
#  @seemethod "getReplicates" with \code{nbrOfReplicates=3}.
# }
#
# \arguments{
#   \item{...}{Arguments accepted by @seemethod "getReplicates".}
# }
#
# \value{
#   Returns a @list or a @matrix.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReplicates".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getTriplicates", "GeneGroups", function(this, ...) {
  getReplicates(this, nbrOfReplicates=3, ...);
}) # getTriplicates()



#########################################################################/**
# @RdocMethod getQuadruplicates
#
# @title "Gets all quadruplicated spots"
#
# @synopsis
#
# \description{
#  @get "title". This is a shortcut for calling
#  @seemethod "getReplicates" with \code{nbrOfReplicates=4}.
# }
#
# \arguments{
#   \item{...}{Arguments accepted by @seemethod "getReplicates".}
# }
#
# \value{
#   Returns a @list or a @matrix.
# }
#
# @author
#
# \seealso{
#   @seemethod "getReplicates".
#   @seeclass
# }
#*/#########################################################################
setMethodS3("getQuadruplicates", "GeneGroups", function(this, ...) {
  getReplicates(this, nbrOfReplicates=4, ...);
}) # getQuadruplicates()



############################################################################
# HISTORY:
# 2004-01-12
# o Added Rdoc comments to getDuplicates(), getTriplicates(), 
#   getQuadruplicates(), getSizes(), nbrOfGroups() and getNames().
# 2003-07-21
# o Added getDuplicates(), getTriplicates(), getQuadruplicates(), and the
#   general getReplicates().
# o Moved into its own source file.
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
