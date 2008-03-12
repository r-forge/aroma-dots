###########################################################################/**
# @RdocClass AffymetrixCelSetTuple
#
# @title "The AffymetrixCelSetTuple class"
#
# \description{
#  @classhierarchy
# }
# 
# @synopsis
#
# \arguments{
#   \item{csList}{A single or @list of @see "AffymetrixCelSet":s.}
#   \item{tags}{A @character @vector of tags.}
#   \item{...}{Not used.}
#   \item{.setClass}{The name of the class of the input set.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
# 
# \seealso{
# }
#*/###########################################################################
setConstructorS3("AffymetrixCelSetTuple", function(csList=NULL, tags="*", ..., .setClass="AffymetrixCelSet") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'csList':
  if (!is.null(csList)) {
    if (!is.list(csList)) {
      csList <- list(csList);
    }

    for (cs in csList) {
      if (!inherits(cs, .setClass)) {
        throw("Argument 'csList' contains a non-", .setClass, ": ", 
                                                            class(cs)[1]);
      }
    }
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }

  extend(Object(), "AffymetrixCelSetTuple",
    .csList = csList,
    .tags = tags
  )
})


setMethodS3("clearCache", "AffymetrixCelSetTuple", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.cesList))
    clearCache(this$.cesList);

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("byPath", "AffymetrixCelSetTuple", function(static, path, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'path':
  path <- Arguments$getReadablePath(path, mustExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Get the corresponding CEL set class
  className <- class(static)[1];
  className <- gsub("Tuple$", "", className);
  clazz <- Class$forName(className);

  verbose && enter(verbose, "Creating ", class(static)[1]);
  verbose && cat(verbose, "Path: ", path);

  # Get all subdirectories
  verbose && enter(verbose, "Scanning path for subdirectories");
  dirs <- list.files(path=path, full.names=TRUE);
  dirs <- dirs[sapply(dirs, FUN=isDirectory)];
  verbose && print(verbose, dirs);
  verbose && exit(verbose);

  # Check which corresponds to chip types, i.e. has CDF files
  verbose && enter(verbose, "Keeping those with names of chip types");
  dirs <- sapply(dirs, FUN=function(dir) {
    chipType <- basename(dir);
    pathname <- AffymetrixCdfFile$findByChipType(chipType);
    if (is.null(pathname))
      dir <- NA;
    dir;
  })
  dirs <- dirs[!is.na(dirs)];
  verbose && print(verbose, dirs);
  verbose && exit(verbose);

  # Define CEL sets for each directory
  verbose && enter(verbose, "Defining list of ", className, ":s");
  csList <- base::lapply(dirs, FUN=function(dir) {
    clazz$fromFiles(dir, ..., verbose=less(verbose));
  })
  names(csList) <- basename(names(csList));
  verbose && str(verbose, dirs);
  verbose && exit(verbose);

  verbose && exit(verbose);

  newInstance(static, csList, ...);
}, static=TRUE)



setMethodS3("as.character", "AffymetrixCelSetTuple", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip types:", paste(getChipTypes(this), collapse=", ")));
  csList <- getListOfSets(this);
  for (kk in seq(along=csList)) {
    cs <- csList[[kk]];
    s <- c(s, as.character(cs));
  }
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)



###########################################################################/**
# @RdocMethod extract
#
# @title "Extracts a subset AffymetrixCelSetTuple"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of arrays to be extracted.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @see "AffymetrixCelSetTuple".
# }
#
# \details{
#   If no matching arrays are available for a given chip type, that chip type
#   is excluded in the returned set tuple.  This guarantees that a set tuple
#   always contains existing arrays.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extract", "AffymetrixCelSetTuple", function(this, arrays, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting subset of arrays");

  # Identify the array indices for each chip type
  arrayTable <- getTableOfArrays(this);
  rr <- indexOfArrays(this, arrays);
  verbose && print(verbose, rr);
  arrayTable <- arrayTable[rr,,drop=FALSE];
  cc <- which(apply(arrayTable, MARGIN=2, FUN=function(idxs) {
    any(!is.na(idxs));
  }));
  verbose && print(verbose, cc);
  arrayTable <- arrayTable[,cc,drop=FALSE];

  verbose && print(verbose, arrayTable);
  
  # Extract only those arrays
  res <- clone(this);
  csList <- getListOfSets(this)[cc];
  verbose && str(verbose, csList);

  for (kk in seq(along=csList)) {
    cs <- csList[[kk]];
    idxs <- arrayTable[,kk];
    idxs <- na.omit(idxs);
    cs <- extract(cs, idxs);
    csList[[kk]] <- cs;
  }
  res$.csList <- csList;

  verbose && str(verbose, csList);

  verbose && exit(verbose);

  res;
})




setMethodS3("getName", "AffymetrixCelSetTuple", function(this, collapse="+", ...) {
  name <- getAlias(this);

  if (is.null(name)) {
    # Get name of chip-effect sets
    csList <- getListOfSets(this);
  
    # Get names
    names <- lapply(csList, FUN=getName);
    names <- unlist(names, use.names=FALSE);
  
    # Merge names
    names <- mergeByCommonTails(names, collapse=collapse);

    name <- names;
  }

  name;
})


setMethodS3("getAlias", "AffymetrixCelSetTuple", function(this, ...) {
  this$.alias;
})


setMethodS3("setAlias", "AffymetrixCelSetTuple", function(this, alias=NULL, ...) {
  # Argument 'alias':
  alias <- Arguments$getCharacter(alias);
  this$.alias <- alias;
  invisible(this);
})


setMethodS3("getAsteriskTags", "AffymetrixCelSetTuple", function(this, ...) {
  "";
}, protected=TRUE)


setMethodS3("setTags", "AffymetrixCelSetTuple", function(this, tags=NULL, ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }

  this$.tags <- tags;

  invisible(this);
})


setMethodS3("getTags", "AffymetrixCelSetTuple", function(this, collapse=NULL, ...) {
  # Get tags of chip-effect set
  csList <- getListOfSets(this);

  # Get data set tags
  tags <- lapply(csList, FUN=getTags);

  # Keep common tags
  tags <- getCommonListElements(tags);
  tags <- tags[[1]];
  tags <- unlist(tags, use.names=FALSE);

  # Add optional tuple tags
  tags <- c(tags, this$.tags);

  # Update asterisk tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Remove empty tags
  tags <- tags[nchar(tags) > 0];

  # Remove duplicated tags 
  tags <- locallyUnique(tags);

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    if (length(tags) > 0)
      tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("getFullName", "AffymetrixCelSetTuple", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



setMethodS3("getListOfSets", "AffymetrixCelSetTuple", function(this, ...) {
  sets <- this$.csList;
  if (is.null(names(sets))) {
    names(sets) <- sapply(sets, FUN=function(set) {
      cdf <- getCdf(set);
      getChipType(cdf, fullname=FALSE);
    });
    this$.csList <- sets;
  }
  sets;
})



###########################################################################/**
# @RdocMethod getChipTypes
#
# @title "Gets a label for all chip types merged"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{merge}{If @TRUE, the chip types are merged into a single string
#     compressed such that only non-common parts are replicated.}
#   \item{collapse}{The string used to merge chip types.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipTypes", "AffymetrixCelSetTuple", function(this, merge=FALSE, collapse="+", ...) {
  cdfList <- getListOfCdfs(this);
  chipTypes <- sapply(cdfList, FUN=getChipType, fullname=FALSE);

  # Invariant for order
#  chipTypes <- sort(chipTypes);

  # Merge to a single string?
  if (merge)
    chipTypes <- mergeByCommonTails(chipTypes, collapse=collapse);

  chipTypes;
})


###########################################################################/**
# @RdocMethod nbrOfChipTypes
#
# @title "Gets the number of chip types"
#
# \description{
#  @get "title" used in the model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfChipTypes", "AffymetrixCelSetTuple", function(this, ...) {
  length(getChipTypes(this, ...));
})


setMethodS3("getListOfCdfs", "AffymetrixCelSetTuple", function(this, ...) {
  csList <- getListOfSets(this);
  lapply(csList, FUN=getCdf);
}, private=TRUE)





###########################################################################/**
# @RdocMethod getTableOfArrays
#
# @title "Gets a table of arrays"
#
# \description{
#  @get "title" showing their availability across chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a \eqn{NxK} @matrix of @integers where \eqn{N} is the total number 
#  of arrays and \eqn{K} is the number of chip types in the model.  The row 
#  names are the names of the arrays, and the column names are the chip types.
#  If data is available for array \eqn{n} and chip type \eqn{k}, cell 
#  \eqn{(n,k)} has value \eqn{n}, otherwise @NA.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTableOfArrays", "AffymetrixCelSetTuple", function(this, ...) {
  csList <- getListOfSets(this);

  # Get all chip types for this data set
  chipTypes <- getChipTypes(this);
  nbrOfChipTypes <- length(csList);

  # Get all sample names
  names <- lapply(csList, FUN=getNames);
  names(names) <- chipTypes;

  # Get all sample names
  allNames <- unlist(names, use.names=FALSE);

  # Get all unique names
  allNames <- unique(allNames);
  lens <- sapply(names, FUN=length);
  nbrOfArrays <- max(lens, length(allNames));
  allNames <- rep(allNames, length.out=nbrOfArrays);

  # Create table of arrays
  nbrOfArrays <- length(allNames);
  X <- matrix(NA, nrow=nbrOfArrays, ncol=nbrOfChipTypes);
  dimnames(X) <- list(allNames, chipTypes);
  for (chipType in chipTypes) {
    namesCT <- names[[chipType]];
    for (rr in seq(length=nrow(X))) {
      name <- rownames(X)[rr];
      idx <- match(name, namesCT);
      if (!is.na(idx)) {
        X[idx,chipType] <- rr;
        namesCT[idx] <- NA;
      }
    }
  }

  X;
}, protected=TRUE)


setMethodS3("getNames", "AffymetrixCelSetTuple", function(this, ...) {
  rownames(getTableOfArrays(this, ...));
})


setMethodS3("getFullNames", "AffymetrixCelSetTuple", function(this, arrays=NULL, exclude=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getFullNameOfTuple <- function(cfList) {
    # Get sample name
    first <- which(!sapply(cfList, FUN=is.null))[1];
    name <- getName(cfList[[first]]);
  
    # Get chip-effect tags *common* across chip types
    tags <- lapply(cfList, FUN=function(ce) {
      if (is.null(ce)) NULL else getTags(ce);
    });
    tags <- base::lapply(tags, setdiff, exclude);
    tags <- getCommonListElements(tags);
    tags <- tags[[1]];
    tags <- unlist(tags, use.names=FALSE);
    tags <- locallyUnique(tags);
    
    fullname <- paste(c(name, tags), collapse=",");
    
    fullname;
  } # getFullNameOfTuple()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq_len(nbrOfArrays(this));
} else {
    arrays <- Arguments$getIndices(arrays, range=c(1, nbrOfArrays(this)));
  }

  # Argument 'exclude':
  exclude <- Arguments$getCharacters(exclude);

  
  fullnames <- c();
  for (kk in arrays) {
    cfList <- getArrayTuple(this, array=kk, ...);
    fullname <- getFullNameOfTuple(cfList);
    fullnames <- c(fullnames, fullname);
  }

  fullnames;
})



###########################################################################/**
# @RdocMethod getArrays
#
# @title "Gets the names of the arrays"
#
# \description{
#  @get "title" available in the model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getArrays", "AffymetrixCelSetTuple", function(this, ...) {
  arrays <- getNames(this, ...);
  names(arrays) <- getFullNames(this);
  arrays;
})



###########################################################################/**
# @RdocMethod indexOfArrays
#
# @title "Gets the indices of the arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @character @vector of arrays names.
#     If @NULL, all arrays are considered.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("indexOfArrays", "AffymetrixCelSetTuple", function(this, arrays=NULL, ...) {
  allNames <- getNames(this);

  # Argument 'arrays':
  if (is.null(arrays)) {
    arrays <- seq(along=allNames);
  } else if (is.numeric(arrays)) {
    arrays <- Arguments$getIndices(arrays, range=c(1,length(allNames)));
  } else {
    missing <- which(!(arrays %in% allNames));
    if (length(missing) > 0) {
      missing <- paste(arrays[missing], collapse=", ");
      throw("Argument 'arrays' contains unknown arrays: ", missing);
    }
    arrays <- match(arrays, allNames);
  }

  arrays;
}, private=TRUE)



###########################################################################/**
# @RdocMethod nbrOfArrays
#
# @title "Gets the number of arrays"
#
# \description{
#  @get "title" used in the model.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @integer.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("nbrOfArrays", "AffymetrixCelSetTuple", function(this, ...) {
  length(getNames(this, ...));
})


###########################################################################/**
# @RdocMethod getArrayTuple
#
# @title "Gets arrays across chip types for one sample"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{array}{Sample of interest.}
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a named @list of @see "AffymetrixCelSet":s.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getArrayTuple", "AffymetrixCelSetTuple", function(this, array, ...) {
  # Get table of arrays
  arrayTable <- getTableOfArrays(this);

  # Argument 'array':
  if (is.numeric(array)) {
    array <- Arguments$getIndex(array, range=c(1, nrow(arrayTable)));
  } else {
    array <- Arguments$getCharacter(array);
    arrayNames <- rownames(arrayTable);
    if (!array %in% arrayNames)
      throw("Argument 'array' refers to a non-existing array: ", array);
  }


  # Get tuple of array indices across chip types
  idxs <- arrayTable[array,];

  # Allocated array tuple
  chipTypes <- colnames(arrayTable);
  tuple <- vector("list", length(chipTypes));

  # Extract tuple
  csList <- getListOfSets(this);
  for (kk in seq(along=csList)) {
    idx <- idxs[kk];
    if (is.na(idx))
      next;
    cs <- csList[[kk]];
    tuple[[kk]] <- getFile(cs, idx);
  }
  names(tuple) <- chipTypes;

  tuple;
})


setMethodS3("asMatrixOfFiles", "AffymetrixCelSetTuple", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting matrix of files");

  csList <- getListOfSets(this);

  # Identify the array indices for each chip type
  arrayTable <- getTableOfArrays(this);
  verbose && print(verbose, arrayTable);

  arrayOfFiles <- rep(NA, length(arrayTable));
  arrayOfFiles <- as.list(arrayOfFiles);
  dim(arrayOfFiles) <- dim(arrayTable);
  dimnames(arrayOfFiles) <- dimnames(arrayTable);

  for (cc in seq(length=ncol(arrayOfFiles))) {
    files <- as.list(csList[[cc]]);
    files <- files[arrayTable[,cc]];
    arrayOfFiles[,cc] <- files;
  }
  rm(files, arrayTable, csList);

  verbose && exit(verbose);

  arrayOfFiles;
}, protected=TRUE)


##############################################################################
# HISTORY:
# 2008-03-11
# o Renamed getTuple() to getArrayTuple().
# 2007-03-29
# o Added asMatrixOfFiles().
# 2007-03-20
# o Now getArrays() returns a named list where the names are the result from
#   getFullNames().
# 2007-03-19
# o TODO: Handle replicated sample names. It is not clear how this should be
#   done.
# o Created from GladModel.R.
##############################################################################
