###########################################################################/**
# @RdocClass CopyNumberSegmentationModel
#
# @title "The CopyNumberSegmentationModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a copy-number segmentation model.
# }
# 
# @synopsis
#
# \arguments{
#   \item{cesTuple}{A @see "ChipEffectSetTuple".}
#   \item{refTuple}{An optional @see "ChipEffectSetTuple" for pairwise 
#     comparisons.}
#   \item{tags}{A @character @vector of tags.}
#   \item{genome}{A @character string specifying what genome is process.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Requirements}{
#   This class requires genome information annotation files for 
#   every chip type.
# }
#
# @author
#*/###########################################################################
setConstructorS3("CopyNumberSegmentationModel", function(cesTuple=NULL, refTuple=NULL, tags="", genome="Human", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cesTuple':
  if (!is.null(cesTuple)) {
    # Coerce to ChipEffectSetTuple
    if (!inherits(cesTuple, "ChipEffectSetTuple")) {
      cesTuple <- ChipEffectSetTuple(cesTuple);
    }

    # Assert that we are dealing with CnChipEffectSet:s
    for (ces in getListOfSets(cesTuple)) {
      # Assert special properties for CnChipEffectSet:s AD HOC /HB 2006-12-20
      if (inherits(ces, "CnChipEffectSet")) {
        # Currently only total copy-number estimates are accepted
        if (!ces$combineAlleles) {
          throw("Unsupported copy-number chip effects. Currently only total copy-number estimates are supported: ces$combineAlleles == FALSE");
        }
      }
    }
  }

  # Argument 'refTuple':
  if (!is.null(refTuple)) {
    # Coerce to ChipEffectSetTuple
    if (!inherits(refTuple, "ChipEffectSetTuple")) {
      refTuple <- ChipEffectSetTuple(refTuple);
    }

    # Assert the same number of chip types in test and reference set
    if (nbrOfChipTypes(refTuple) != nbrOfChipTypes(cesTuple)) {
      throw("The number of reference sets does not match the number of test sets: ", nbrOfChipTypes(refTuple), " != ", nbrOfChipTypes(cesTuple));
    }

    # Assert that we are dealing with CnChipEffectSet:s
    for (ref in getListOfSets(refTuple)) {
      # Assert special properties for CnChipEffectSet:s AD HOC /HB 2006-12-20
      if (inherits(ref, "CnChipEffectSet")) {
        # Currently only total copy-number estimates are accepted
        if (!ref$combineAlleles) {
          throw("Unsupported copy-number chip effects. Currently only total copy-number estimates are supported: ref$combineAlleles == FALSE");
        }
      }
    }

    # Validate consistency between the chip-effect sets and the reference files
    cesList <- getListOfSets(cesTuple);
    refList <- getListOfSets(refTuple);
    for (kk in seq(along=cesList)) {
      ces <- cesList[[kk]];
      ref <- refList[[kk]];

      # Assert that the reference and the test sets are of the same size
      if (nbrOfArrays(ref) != nbrOfArrays(ces)) {
        throw("The number of reference files does not match the number of test files: ", nbrOfArrays(ref), " != ", nbrOfArrays(ces));
      }

      # Assert special properties for CnChipEffectSet:s AD HOC /HB 2006-12-20
      if (inherits(ces, "CnChipEffectSet")) {
        if (ref$combineAlleles != ces$combineAlleles) {
           throw("The reference chip effects are not compatible with the chip-effect set. One is combining the alleles the other is not.");
        }

        if (ref$mergeStrands != ces$mergeStrands) {
           throw("The reference chip effects are not compatible with the chip-effect set. One is merging the strands the other is not.");
        }
      }

      # Assert that the reference files are compatible with the test files
      for (jj in seq(length=nbrOfArrays(ces))) {
        cf <- getFile(ces, jj);
        rf <- getFile(ref, jj);
        if (!inherits(rf, class(cf)[1])) {
          throw(class(ref)[1], " #", kk, " of argument 'refTuple' contains file #", jj, ", that is not of the same class as the paired test file: ", class(rf)[1], " !inherits from ", class(cf)[1]);
        }
      }
    } # for (kk in ...)
  }

  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
  }


  this <- extend(Object(), "CopyNumberSegmentationModel",
    .cesTuple = cesTuple,
    .refTuple = refTuple,
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome
  );

  # Validate?
  if (!is.null(this$.cesTuple)) {
    # Validate genome
    pathname <- getGenomeFile(this);
  }

  this;
})


setMethodS3("as.character", "CopyNumberSegmentationModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", paste(getTags(this), collapse=",")));
  s <- c(s, paste("Chip type (virtual):", getChipType(this)));
  s <- c(s, sprintf("Path: %s", getPath(this)));
  nbrOfChipTypes <- nbrOfChipTypes(this);
  s <- c(s, sprintf("Number of chip types: %d", nbrOfChipTypes));
  s <- c(s, "Chip-effect set & reference file pairs:");
  cesList <- getListOfSets(getSetTuple(this));
  refList <- getRefSetTuple(this);
  if (!is.null(refList))
    refList <- getListOfSets(refList);
  for (kk in seq(along=cesList)) {
    s <- c(s, sprintf("Chip type #%d of %d ('%s'):", kk, nbrOfChipTypes, names(cesList)[kk]));
    s <- c(s, "Chip-effect set:");
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    s <- c(s, as.character(ces));
    s <- c(s, "Reference file:");
    if (is.null(ref)) {
      s <- c(s, "<average across arrays>");
    } else {
      s <- c(s, as.character(ref));
    }
  }
#  s <- c(s, "Genome information:", as.character(getGenomeInformation(this)));
  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("clearCache", "CopyNumberSegmentationModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.cesTuple))
    clearCache(this$.cesTuple);
  if (!is.null(this$.refTuple))
    clearCache(this$.refTuple);

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("getSetTuple", "CopyNumberSegmentationModel", function(this, ...) {
  this$.cesTuple;
}, protected=TRUE)

setMethodS3("getRefSetTuple", "CopyNumberSegmentationModel", function(this, ...) {
  this$.refTuple;
}, protected=TRUE)


setMethodS3("getListOfChipEffectSets", "CopyNumberSegmentationModel", function(this, ...) {
  cesTuple <- getSetTuple(this);
  getListOfSets(cesTuple);
})


setMethodS3("getListOfChipEffects", "CopyNumberSegmentationModel", function(this, ...) {
  getListOfChipEffectSets(this, ...);
}, private=TRUE, deprecated=TRUE)



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
setMethodS3("nbrOfChipTypes", "CopyNumberSegmentationModel", function(this, ...) {
  nbrOfChipTypes(getSetTuple(this), ...);
})



setMethodS3("getListOfCdfs", "CopyNumberSegmentationModel", function(this, ...) {
  getListOfCdfs(getSetTuple(this), ...);
}, private=TRUE)


setMethodS3("getChipTypes", "CopyNumberSegmentationModel", function(this, ...) {
  getChipTypes(getSetTuple(this), ...);
})


###########################################################################/**
# @RdocMethod getChipType
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
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipType", "CopyNumberSegmentationModel", function(this, ...) {
  getChipTypes(this, merge=TRUE, ...);
})


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
setMethodS3("getTableOfArrays", "CopyNumberSegmentationModel", function(this, ...) {
  getTableOfArrays(getSetTuple(this), ...);
})


setMethodS3("getNames", "CopyNumberSegmentationModel", function(this, ...) {
  rownames(getTableOfArrays(this, ...));
})


setMethodS3("getFullNames", "CopyNumberSegmentationModel", function(this, ...) {
  getFullNames(getSetTuple(this), ...);
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
setMethodS3("getArrays", "CopyNumberSegmentationModel", function(this, ...) {
  getNames(this, ...);
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
setMethodS3("indexOfArrays", "CopyNumberSegmentationModel", function(this, arrays=NULL, ...) {
  indexOfArrays(getSetTuple(this), arrays=arrays, ...);
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
setMethodS3("nbrOfArrays", "CopyNumberSegmentationModel", function(this, ...) {
  length(getNames(this, ...));
})


setMethodS3("getName", "CopyNumberSegmentationModel", function(this, collapse="+", ...) {
  name <- getAlias(this);

  if (is.null(name))
    name <- getName(getSetTuple(this), ...);

  name;
})

setMethodS3("getAlias", "CopyNumberSegmentationModel", function(this, ...) {
  this$.alias;
})


setMethodS3("setAlias", "CopyNumberSegmentationModel", function(this, alias=NULL, ...) {
  # Argument 'alias':
  alias <- Arguments$getCharacter(alias);
  this$.alias <- alias;
  invisible(this);
})





setMethodS3("getAsteriskTag", "CopyNumberSegmentationModel", function(this, ...) {
  # Default '*' tag is the abbreviation from upper-case letters only,
  # e.g. "FooHooMooModel" gives "FHM". 
  tag <- class(this)[1];
  tag <- gsub("Model$", "", tag);
  tag <- strsplit(tag, split="")[[1]];
  tagUC <- toupper(tag);
  keep <- (tag == tagUC);
  tag <- tag[keep];
  tag <- paste(tag, collapse=""); 

  # Paired?
  refTuple <- getRefSetTuple(this);
  if (!is.null(refTuple)) {
#    tag <- c(tag, "paired");
#    tag <- paste(tag, collapse=",");
  }

  tag;
}, protected=TRUE)



setMethodS3("getTags", "CopyNumberSegmentationModel", function(this, collapse=NULL, ...) {
  tags <- getTags(getSetTuple(this), ...);

  # Add model tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTag(this);

  # Keep non-empty tags
  tags <- tags[nchar(tags) > 0];

  # Get unique tags
  tags <- unique(tags);

  tags <- paste(tags, collapse=collapse);
  if (length(tags) == 0)
    tags <- NULL;

  tags;
})


setMethodS3("getFullName", "CopyNumberSegmentationModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "CopyNumberSegmentationModel", function(this, ...) {
  # Example: gladData/.
  sprintf("%sData", tolower(getAsteriskTag(this)));
}, private=TRUE)


setMethodS3("getPath", "CopyNumberSegmentationModel", function(this, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Chip type
  chipType <- getChipType(this);

	  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})

###########################################################################/**
# @RdocMethod getChromosomes
#
# @title "Gets the chromosomes available"
#
# \description{
#  @get "title".
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
setMethodS3("getChromosomes", "CopyNumberSegmentationModel", function(this, ...) {
  gis <- getListOfGenomeInformations(this);
  chromosomes <- lapply(gis, getChromosomes);
  chromosomes <- unlist(chromosomes, use.names=TRUE);
  chromosomes <- sort(unique(chromosomes));
  chromosomes;
})


setMethodS3("getListOfGenomeInformations", "CopyNumberSegmentationModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving genome informations");
  cdfList <- getListOfCdfs(this);
  giList <- lapply(cdfList, getGenomeInformation, verbose=less(verbose));
  verbose && exit(verbose);

  giList;
})


setMethodS3("getChipEffectFiles", "CopyNumberSegmentationModel", function(this, ...) {
  setTuple <- getSetTuple(this);
  getTuple(setTuple, ...);
})



setMethodS3("getListOfReferences", "CopyNumberSegmentationModel", function(this, ..., .calculateAvgs=TRUE) {
  refList <- getRefSetTuple(this);
  if (!is.null(refList)) {
    refList <- getListOfSets(refList);
  } else {
    refList <- vector("list", nbrOfChipTypes(this));
  }
  refList;
})


setMethodS3("getReferenceName", "CopyNumberSegmentationModel", function(this, collapse="+", ...) {
  # Get name of the reference files
  refList <- getListOfReferences(this);

  # Get names
  names <- lapply(refList, FUN=function(file) {
    if (is.null(file)) NULL else getName(file);
  });
  names <- unlist(names, use.names=FALSE);

  # Merge names
  names <- mergeByCommonTails(names, collapse=collapse);

  names;
}, private=TRUE)



setMethodS3("getReferenceFiles", "CopyNumberSegmentationModel", function(this, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Retrieving reference data files");

  cesList <- getListOfChipEffects(this);
  refList <- getListOfReferences(this);
  for (kk in seq(length=nbrOfChipTypes(this))) {
    ces <- cesList[[kk]];
    ref <- refList[[kk]];
    if (force || is.null(ref)) {
      if (force) {
        verbose && cat(verbose, "Forced recalculation requested.");
      } else {
        verbose && cat(verbose, "No reference available.");
      }
      verbose && enter(verbose, "Calculating average chip effects");
      ref <- getAverageFile(ces, force=force, verbose=less(verbose));
      verbose && exit(verbose);
      refList[[kk]] <- ref;
    }
  }
  verbose && exit(verbose);

  refList;
})


setMethodS3("getGenome", "CopyNumberSegmentationModel", function(this, ...) {
  this$.genome;
})


setMethodS3("getGenomeFile", "CopyNumberSegmentationModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  fullname <- getGenome(this);
  pattern <- sprintf("^%s.*,chromosomes.txt$", fullname);

  # 1. Search in the regular places
  pathname <- findAnnotationData(name=fullname, set="genomes", 
                            pattern=pattern, ..., verbose=less(verbose, 10));

  # 2. As a backup, search in the <pkg>/annotationData/ directory
  if (is.null(pathname)) {
    verbose && enter(verbose, "Search among package's annotationData/");
    path <- system.file("annotationData", package="aroma.affymetrix");
    verbose && cat(verbose, "Path: ", path);
    pathname <- findAnnotationData(name=fullname, set="genomes", 
                pattern=pattern, ..., paths=path, verbose=less(verbose, 10));
    verbose && exit(verbose);
  }

  if (is.null(pathname)) {
    throw("Failed to locate a genome annotation data file: ", fullname);
  }

  pathname;
}, protected=TRUE)


setMethodS3("setGenome", "CopyNumberSegmentationModel", function(this, genome, tags=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  oldGenome <- this$.genome;

  fullname <- paste(c(genome, tags), collapse=",");
  verbose && cat(verbose, "Fullname: ", fullname);

  # Verify that there is an existing genome file
  tryCatch({
    this$.genome <- fullname;
    pathname <- getGenomeFile(this, verbose=less(verbose, 10));
  }, error = function(ex) {
    this$.genome <- oldGenome;
    throw(ex$message);
  })

  invisible(oldGenome);
})



setMethodS3("getGenomeData", "CopyNumberSegmentationModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading genome chromosome annotation file");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get genome annotation data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for the file");
  # Search annotationData/genomes/
  pathname <- getGenomeFile(this, verbose=less(verbose, 10));
  verbose && exit(verbose);

  verbose && enter(verbose, "Reading data file");
  verbose && cat(verbose, "Pathname: ", pathname);
  data <- readTable(pathname, header=TRUE, 
                            colClasses=c(nbrOfBases="integer"), row.names=1);
  verbose && exit(verbose);

  verbose && enter(verbose, "Translating chromosome names");
  chromosomes <- row.names(data);
  map <- c("X"=23, "Y"=24, "Z"=25);
  for (kk in seq(along=map)) {
    chromosomes <- gsub(names(map)[kk], map[kk], chromosomes, fixed=TRUE);
  }
  row.names(data) <- chromosomes;
  verbose && exit(verbose);

  verbose && exit(verbose);

  data;
}, protected=TRUE)


setMethodS3("getRawCnData", "CopyNumberSegmentationModel", function(this, ceList, refList, chromosome, units=NULL, reorder=TRUE, ..., maxNAFraction=1/8, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ceList':
  if (!inherits(ceList, "list")) {
    throw("Argument 'ceList' is not a list: ", class(ceList)[1]);
  } else {
    for (kk in seq(along=ceList)) {
      if (!inherits(ceList[[kk]], "ChipEffectFile")) {
        throw("Argument 'ceList' contains a non-ChipEffectFile: ", 
                                                     class(ceList[[kk]])[1]);
      }
    }
  }

  # Argument 'refList':
  if (!inherits(refList, "list")) {
    throw("Argument 'refList' is not a list: ", class(refList)[1]);
  } else {
    for (kk in seq(along=refList)) {
      if (!inherits(refList[[kk]], "ChipEffectFile")) {
        throw("Argument 'refList' contains a non-ChipEffectFile: ", 
                                                    class(refList[[kk]])[1]);
      }
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Retrieving raw CN data");

  # Data set attributes
  chipTypes <- getChipTypes(this);
  arrayNames <- getArrays(this);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (x, M, stddev, chiptype, unit) from all chip types
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving relative chip-effect estimates");
  # Get the chip types as a factor
  chipTypes <- as.factor(chipTypes);
  df <- NULL;
  for (kk in seq(along=chipTypes)) {
    chipType <- chipTypes[kk];
    verbose && enter(verbose, "Chip type: ", chipType);
    ce <- ceList[[kk]];
    if (!is.null(ce)) {
      ref <- refList[[kk]];
      df0 <- getXAM(ce, other=ref, chromosome=chromosome, units=units, verbose=less(verbose));
      df0 <- df0[,c("x", "M")];
      verbose && cat(verbose, "Number of units: ", nrow(df0));

      # Estimate the std dev of the raw log2(CN).  [only if ref is average across arrays]
      units0 <- as.integer(rownames(df0));
      # Get (mu, sigma) of theta (estimated across all arrays).
      data <- getDataFlat(ref, units=units0, verbose=less(verbose));
      # Number of arrays (for each unit)
      n <- readCel(getPathname(ref), indices=data$cell, readIntensities=FALSE, readPixels=TRUE)$pixels;

      # Use Gauss' approximation (since mu and sigma are on the 
      # intensity scale)
      sdM <- log2(exp(1)) * sqrt(1+1/n) * data$sdTheta / data$theta;
      rm(n);

      verbose && enter(verbose, "Scanning for non-finite values");
      n <- sum(!is.finite(df0[,"M"]));
      fraction <- n / nrow(df0);
      verbose && printf(verbose, "Number of non-finite values: %d (%.1f%%)\n", 
                                                             n, 100*fraction);

      # Sanity check
      if (fraction > maxNAFraction) {
        throw(sprintf("Something is wrong with the data. Too many non-finite values: %d (%.1f%% > %.1f%%)", as.integer(n), 100*fraction, 100*maxNAFraction));
      }
      verbose && exit(verbose);
  
      # Append SD, chip type, and CDF units.
      df0 <- cbind(df0, sdTheta=data$sdTheta, sdM=sdM, chipType=rep(chipType, length=length(units0)), unit=units0);
      rm(data);
  
      df <- rbind(df, df0);
      colnames(df) <- colnames(df0);
      rm(df0, units0);
    } else {
      verbose && cat(verbose, "No chip-effect estimates available: ", arrayNames[kk]);
    }

    # Garbage collect
    gc <- gc();
    
    verbose && exit(verbose);
  } # for (kk in ...)
  verbose && exit(verbose);
  

  if (reorder) {
    verbose && enter(verbose, "Re-order by physical position");
    df <- df[order(df[,"x"]),];
    rownames(df) <- NULL;
    nbrOfUnits <- nrow(df);
    verbose && exit(verbose);
  }

  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);

  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nbrOfUnits));
  verbose && exit(verbose);

  df;
}, protected=TRUE);




###########################################################################/**
# @RdocMethod fitOne
#
# @title "Fits the copy-number segmentation model for one chromosome in one sample"
#
# \description{
#  @get "title".
#
#  \emph{This is an abstract method that has to be implemented in a subclass.}
# }
#
# @synopsis
#
# \arguments{
#   \item{data}{A @data.frame with columns \code{M} (log-ratio) and 
#      \code{x} (locus position).
#   }
#   \item{chromosome}{An @integer specifying the index of the chromosome to
#      be fitted.}
#   \item{...}{Additional arguments passed down to the internal fit function.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an object returned by internal fit function.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitOne", "CopyNumberSegmentationModel", abstract=TRUE);



###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{arrays}{A @vector of array indices specifying which arrays to
#    be considered.  If @NULL, all are processed.}
#   \item{chromosome}{A @vector of chromosomes indices specifying which
#     chromosomes to be considered.  If @NULL, all are processed.}
#   \item{force}{If @FALSE, the model will not be fitted again if it was
#     already fitted.}
#   \item{...}{Not used.}
#   \item{.retResults}{If @TRUE, CBS fit structures are returned for each
#     fitted array and chromosome.}
#   \item{verbose}{A @logical or @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a named @list of named @lists.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fit", "CopyNumberSegmentationModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), force=FALSE, ..., .retResults=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  if (identical(arrays, "fitted")) {
  } else {
    arrays <- indexOfArrays(this, arrays=arrays);
  }

  allChromosomes <- getChromosomes(this);

  # Argument 'chromosomes':
  if (identical(chromosomes, "fitted")) {
  } else if (is.null(chromosomes)) {
    chromosomes <- getChromosomes(this);
  } else if (is.numeric(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes,
                                                range=range(allChromosomes));
##    chromosomes <- as.character(chromosomes);  ## TODO
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, allChromosomes);
  } else if (is.character(chromosomes)) {
    chromosomes <- Arguments$getChromosomes(chromosomes, 
                                                range=range(allChromosomes));
##    chromosomes[chromosomes == "23"] <- "X";   ## TODO
    chromosomes <- intersect(chromosomes, getChromosomes(this));
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- getPath(this);
  mkdirs(path);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving reference chip effects
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Note: Do *not* pass done 'force' to getReferenceFiles(), because then
  # it will calculate the average regardless of reference. /HB 2007-03-24
  refList <- getReferenceFiles(this, verbose=verbose);
  verbose && cat(verbose, "Using references:");
  verbose && print(verbose, refList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get reference annotations
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get reference tags across chip types
  tags <- lapply(refList, getTags);
  # HB 2007-02-19 To fix: Should the name and the tags of average files
  # be replaced?!? We do not get the same names if we specify the average
  # files explicitly or indirectly.
  #  tags <- getCommonListElements(tags);
  tags <- unlist(tags, use.names=FALSE);
  tags <- setdiff(tags, "chipEffects");
  refTags <- tags;
#  verbose && cat(verbose, "Reference tags: ", paste(refTags, collapse=","));

  # Add combined reference name
  refName <- getReferenceName(this);
  refTags <- c(refName, refTags);
  verbose && cat(verbose, "Reference tags: ", paste(refTags, collapse=","));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Chromosome by chromosome
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- list();
  arrayNames <- getNames(this)[arrays];
  nbrOfArrays <- length(arrayNames);
  for (aa in seq(length=nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    ceList <- getChipEffectFiles(this, array=array);
    rfList <- refList;
##    rfList <- lapply(refList, FUN=function(ces) { 
##      if (is.null(ces)) {
##        NULL;
##      } else {
##        getFile(ces, array);
##      }
##    });

    # Get chip-effect tags *common* across chip types
    tags <- lapply(ceList, FUN=function(ce) {
      if (is.null(ce)) NULL else getTags(ce);
    });
    tags <- getCommonListElements(tags);
    tags <- unlist(tags, use.names=FALSE);
    tags <- setdiff(tags, "chipEffects");
    ceTags <- tags;
    verbose && cat(verbose, "Chip-effect tags: ", paste(ceTags, collapse=","));

    res[[arrayName]] <- list();
    for (chr in chromosomes) {
##      chrIdx <- match(chr, c(1:22, "X", "Y"));
      verbose && enter(verbose, 
                          sprintf("Array #%d ('%s') of %d on chromosome %s", 
                                           aa, arrayName, nbrOfArrays, chr));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get pathname
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Add tags chrNN,<reference tags>
      tags <- c(ceTags, sprintf("chr%02d", chr));
      tags <- c(tags, refTags);
      fullname <- paste(c(arrayName, tags), collapse=",");
      filename <- sprintf("%s.xdr", fullname);
      pathname <- filePath(path, filename);

      # Already done?
      if (!force && isFile(pathname)) {
        verbose && enter(verbose, "Loading results from file");
        verbose && cat(verbose, "Pathname: ", pathname);
        fit <- loadObject(pathname);
        verbose && cat(verbose, "Fit object: ", class(fit)[1]);
        verbose && exit(verbose);
      } else {
        # Time the fitting.
        startTime <- processTime();

        timers <- list(total=0, read=0, fit=0, write=0, gc=0);

        tTotal <- processTime();

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Get (x, M, stddev, chiptype, unit) data from all chip types
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tRead <- processTime();
print(list(ceList, rfList));
        data <- getRawCnData(this, ceList=ceList, refList=rfList, 
                     chromosome=chr, ..., force=force, verbose=less(verbose));
        timers$read <- timers$read + (processTime() - tRead);


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Fit segmentation model
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tFit <- processTime();
        fit <- fitOne(this, data=data, chromosome=chr, ..., verbose=less(verbose));
        timers$fit <- timers$fit + (processTime() - tFit);
        rm(data); # Not needed anymore

        verbose && cat(verbose, "Class of fitted object: ", class(fit)[1]);
        verbose && printf(verbose, "Time to fit segmentation model: %.2fmin\n", timers$fit[3]/60);

        verbose && enter(verbose, "Validate that it can be coerced");
        rawCns <- extractRawCopyNumbers(fit);
        nbrOfLoci <- nbrOfLoci(rawCns);
        verbose && print(verbose, rawCns);
        cnRegions <- extractCopyNumberRegions(fit);
        verbose && print(verbose, cnRegions);
        rm(rawCns, cnRegions);  # Not needed anymore
        verbose && exit(verbose);


        # Garbage collection
        tGc <- processTime();
        gc <- gc();
        timers$gc <- timers$gc + (processTime() - tGc);
        verbose && print(verbose, gc);

        verbose && enter(verbose, "Saving to file");
        verbose && cat(verbose, "Pathname: ", pathname);
        tWrite <- processTime();
        saveObject(fit, file=pathname);
        timers$write <- timers$write + (processTime() - tWrite);
        verbose && exit(verbose);

        timers$total <- timers$total + (processTime() - tTotal);

        # Report time profiling
        totalTime <- processTime() - startTime;
  
        if (verbose) {
          t <- totalTime[3];
          printf(verbose, "Total time for chromosome %d: %.2fs == %.2fmin\n", chr, t, t/60);
          t <- totalTime[3]/nbrOfLoci;
          printf(verbose, "Total time per 1000 locus (with %d loci): %.2fs/locus\n", nbrOfLoci, 1000*t);
          # Get distribution of what is spend where
          t <- base::lapply(timers, FUN=function(timer) unname(timer[3]));
          t <- unlist(t);
          t <- 100 * t / t["total"];
          printf(verbose, "Fraction of time spent on different tasks: Fitting: %.1f%%, Reading: %.1f%%, Writing: %.1f%%, Explicit garbage collection: %.1f%%\n", t["fit"], t["read"], t["write"], t["gc"]);
        }
      }

      hookName <- "onFit.CopyNumberSegmentationModel";
      verbose && enter(verbose, sprintf("Calling %s() hooks", hookName));
      callHooks(hookName, fit=fit, chromosome=chr, fullname=fullname);
      verbose && exit(verbose);

      if (.retResults)
        res[[arrayName]][[chr]] <- fit;

      rm(fit);
   
      verbose && exit(verbose);
    } # for (chr in ...)
  } # for (aa in ...)

  invisible(res);
})


setMethodS3("getSetTag", "CopyNumberSegmentationModel", function(this, ...) {
  tolower(getAsteriskTag(this));
}, private=TRUE)

setMethodS3("getReportPath", "CopyNumberSegmentationModel", function(this, ...) {
  rootPath <- "reports";

  # Data set name
  name <- getName(this);

  # Data set tags
  tags <- getTags(this, collapse=",");

  # Get chip type
  chipType <- getChipType(this);

  # Image set
  set <- getSetTag(this);

  # The report path
  path <- filePath(rootPath, name, tags, chipType, set, expandLinks="any");

  path;
}, protected=TRUE)



setMethodS3("getLog2Ratios", "CopyNumberSegmentationModel", function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "getLog2Ratios()");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Obtaining fits (or fit if missing)");
  suppressWarnings({
    res <- fit(this, ..., .retResults=TRUE, verbose=less(verbose,10));
  })
  verbose && exit(verbose);

  verbose && enter(verbose, "Extracting regions from all fits");
  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    # For each chromosome
    for (kk in seq(along=arrayFits)) {
      fit <- arrayFits[[kk]];
      if (!is.null(fit)) {
        verbose && enter(verbose, "Extracting regions for chromosome #", kk);
        suppressWarnings({
          df0 <- getRegions(fit, ...);
        })
        df <- rbind(df, df0);
        verbose && exit(verbose);
      }
    }
    rownames(df) <- seq(length=nrow(df));
    df;
  })
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
}, private=TRUE) # getLog2Ratios()



setMethodS3("getRegions", "CopyNumberSegmentationModel", function(this, ..., url="ucsc", organism="Human", hgVersion="hg17", margin=10, flat=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'url':
  if (identical(url, "ucsc")) {
    # The UCSC browser accepts chromsomes either 'X' or 23.  In other words,
    # we can stick with integers to be more general.
    url <- function(chromosome, start, stop) {
      suppressWarnings({
        start <- as.integer(start);
        stop <- as.integer(stop);
      })
      uri <- "http://genome.ucsc.edu/cgi-bin/hgTracks?clade=vertebrate";
      sprintf("%s&org=%s&db=%s&position=chr%s%%3A%d-%d", uri, organism, hgVersion, chromosome, start, stop);
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extracting regions from all fits");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract the regions for each of the CN model fits (per array)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Obtaining CN model fits (or fit if missing)");
  suppressWarnings({
    res <- fit(this, ..., .retResults=TRUE, verbose=less(verbose,10));
  })
  verbose && exit(verbose);

  res <- lapply(res, FUN=function(arrayFits) {
    df <- NULL;
    # For each chromosome
    for (kk in seq(along=arrayFits)) {
      fit <- arrayFits[[kk]];
      if (!is.null(fit)) {
        verbose && enter(verbose, "Extracting regions for chromosome #", kk);
        suppressWarnings({
#          df0 <- getRegions(fit, ...);
          cnr <- extractCopyNumberRegions(fit, ...);
          df0 <- as.data.frame(cnr);
        })
        df <- rbind(df, df0);
        verbose && exit(verbose);
      }
    }
    rownames(df) <- seq(length=nrow(df));

    verbose && cat(verbose, "Extracted regions:");
    verbose && str(verbose, df);

    # Add URL?
    if (!is.null(url)) {
      chromosome <- df[,"chromosome"];
      start <- df[,"start"];
      stop <- df[,"stop"];
      m <- margin*abs(stop-start);
      start <- start-m;
      start[start < 0] <- 0;
      stop <- stop + m;
      urls <- character(nrow(df));
      for (rr in seq(along=urls)) { 
        urls[rr] <- url(chromosome[rr], start[rr], stop[rr]);
      }
      df <- cbind(df, url=urls);
    }

    df;
  })
  verbose && exit(verbose);

  if (flat) {
    df <- NULL;
    for (kk in seq(along=res)) {
      df <- rbind(df, cbind(sample=names(res)[kk], res[[kk]]));
      res[[kk]] <- NA;
    }
    row.names(df) <- seq(length=nrow(df));
    res <- df;
  }  

  res;
})


setMethodS3("writeRegions", "CopyNumberSegmentationModel", function(this, arrays=NULL, format=c("xls", "wig"), digits=3, ..., oneFile=TRUE, skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'arrays':
  arrays <- indexOfArrays(this, arrays=arrays);

  # Argument 'format':
  format <- match.arg(format);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Setup
  fullname <- getFullName(this);
  arrayNames <- getNames(this);

  path <- getPath(this);
  mkdirs(path);

  if (oneFile) {
    filename <- sprintf("%s,regions.%s", fullname, format); 
    pathname <- filePath(path, filename);
    pathname <- Arguments$getWritablePathname(pathname);
    if (!skip && isFile(pathname)) {
      file.remove(pathname);
    }
  }

  res <- list();
  for (aa in seq(along=arrays)) {
    array <- arrays[aa];
    name <- arrayNames[array];
    verbose && enter(verbose, sprintf("Array #%d ('%s') of %d", 
                                               aa, name, length(arrays)));
    df <- getRegions(this, arrays=array, ..., verbose=less(verbose))[[1]];
    names(df) <- gsub("Smoothing", "log2", names(df));

    if (nrow(df) > 0) {
      if (identical(format, "xls")) {
        # Append column with sample names
        df <- cbind(sample=name, df);
      } else if (identical(format, "wig")) {
        # Write a four column WIG/BED table
        df <- df[,c("chromosome", "start", "stop", "log2")];
  
        # In the UCSC Genome Browser, the maximum length of one element
        # is 10,000,000 bases.  Chop up long regions in shorter contigs.
        verbose && enter(verbose, sprintf("Chopping up too long segment"));
        MAX.LENGTH = 10e6-1;
        start <- df[,"start"];
        stop <- df[,"stop"];
        len <- stop-start;
        tooLong <- which(len > MAX.LENGTH);
        if (length(tooLong) > 0) {
          dfXtra <- NULL;
          for (rr in tooLong) {
            x0 <- start[rr];
            while (x0 < stop[rr]) {
              x1 <- min(x0 + MAX.LENGTH, stop[rr]);
              df1 <- df[rr,];
              df1[,"start"] <- x0;
              df1[,"stop"] <- x1;
              dfXtra <- rbind(dfXtra, df1);          
              x0 <- x1+1;
            }
          }
          df <- df[-tooLong,];
          df <- rbind(df, dfXtra);
          rm(dfXtra);
          row.names(df) <- seq(length=nrow(df));
        }
        verbose && exit(verbose);
        # Make sure the items are ordered correctly
        chrIdx <- as.integer(df[,"chromosome"]);
        o <- order(chrIdx, df[,"start"]);
        df <- df[o,];
  
        # All chromosomes should have prefix 'chr'.
        chrIdx <- as.integer(df[,"chromosome"]);
        ## df[chrIdx == 23,"chromosome"] <- "X"; ## REMOVED 2007-03-15
        df[,"chromosome"] <- paste("chr", df[,"chromosome"], sep="");
      }
  
      # Apply digits
      for (cc in seq(length=ncol(df))) {
        value <- df[,cc];
        if (is.double(value)) {
          df[,cc] <- round(value, digits=digits);
        }
      }
    } # if (nrow(df) > 0)

    if (!oneFile) {
      savename <- name;
      filename <- sprintf("%s,regions.%s", savename, format); 
      pathname <- filePath(path, filename);
      if (!oneFile && !skip && isFile(pathname))
        file.remove(pathname);
    }

    # Writing to file
    verbose && cat(verbose, "Pathname: ", pathname);
    if (identical(format, "xls")) {
      col.names <- (array == arrays[1]);
      suppressWarnings({
        write.table(df, file=pathname, sep="\t", col.names=col.names, row.names=FALSE, quote=FALSE, append=oneFile);
      })
    } else if (identical(format, "wig")) {
      # Write track control
      trackAttr <- c(type="wiggle_0");
      trackAttr <- c(trackAttr, name=sprintf("\"%s\"", name));
      trackAttr <- c(trackAttr, 
                     group=sprintf("\"%s regions\"", getAsteriskTag(this)));
      trackAttr <- c(trackAttr, priority=array);
      trackAttr <- c(trackAttr, graphType="bar");
      trackAttr <- c(trackAttr, visibility="full");
      trackAttr <- c(trackAttr, maxHeightPixels="128:96:64");
      trackAttr <- c(trackAttr, yLineOnOff="on");
# HARD WIRED FOR NOW.  TO DO /hb 2006-11-27
col <- c("117,112,179", "231,41,138");
ylim <- c(-1,1);
      if (!is.null(col)) {
        trackAttr <- c(trackAttr, color=col[1], altColor=col[2]);
      }
      if (!is.null(ylim)) {
        trackAttr <- c(trackAttr, autoScale="off", 
              viewLimits=sprintf("%.2f:%.2f ", ylim[1], ylim[2]));
      }
      trackAttr <- paste(names(trackAttr), trackAttr, sep="=");
      trackAttr <- paste(trackAttr, collapse=" ");
      trackAttr <- paste("track ", trackAttr, "\n", sep="");
      verbose && str(verbose, trackAttr);
      cat(file=pathname, trackAttr, append=oneFile);

      # Write data
      verbose && str(verbose, df);
      write.table(df, file=pathname, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE, append=oneFile);
    }
    verbose && exit(verbose);
    res[[array]] <- df;
  }

  invisible(pathname);
})


setMethodS3("calculateChromosomeStatistics", "CopyNumberSegmentationModel", function(this, arrays=NULL, chromosomes=getChromosomes(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allChromosomes <- getChromosomes(this);

  # Argument 'arrays':
  arrays <- indexOfArrays(this, arrays=arrays);

  # Argument 'chromosomes':
  if (is.null(chromosomes)) {
    chromosomes <- allChromosomes;
  } else {
    unknown <- chromosomes[!(chromosomes %in% allChromosomes)];
    if (length(unknown) > 0) {
      throw("Argument 'chromosomes' contains unknown values: ", paste(unknown, collapse=", "));
    }
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  verbose && enter(verbose, "Estimating number of copies per chromosome");
  nbrOfChromosomes <- length(chromosomes);
  verbose && printf(verbose, "Chromosomes (%d): %s", nbrOfChromosomes, paste(chromosomes, collapse=", "));

  refList <- getReferenceFiles(this, verbose=verbose);
  verbose && cat(verbose, "Using references:");
  verbose && print(verbose, refList);

  res <- list();
  arrayNames <- getNames(this)[arrays];
  nbrOfArrays <- length(arrayNames);

  for (aa in seq(length=nbrOfArrays)) {
    array <- arrays[aa];
    arrayName <- arrayNames[aa];

    ceList <- getChipEffectFiles(this, array=array);

    res[[arrayName]] <- list();
    for (chr in chromosomes) {
      verbose && enter(verbose, 
                          sprintf("Array #%d ('%s') of %d on chromosome %s", 
                                           aa, arrayName, nbrOfArrays, chr));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Get (x, M, stddev, chiptype, unit) data from all chip types
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      data <- getRawCnData(this, ceList=ceList, refList=refList, chromosome=chr, ..., verbose=less(verbose));
      M <- data[,"M"];
      rm(data);

      fit <- list(
        mean = mean(M, na.rm=TRUE),
        sd = sd(M, na.rm=TRUE),
        median = median(M, na.rm=TRUE),
        mad = mad(M, na.rm=TRUE),
        quantiles = quantile(M, probs=seq(0,1,by=0.01), na.rm=TRUE),
        nbrOfLoci = length(M),
        nbrOfNAs = sum(is.na(M))
      );

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Estimate
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Garbage collection
      gc <- gc();
      verbose && print(verbose, gc);

      res[[arrayName]][[chr]] <- fit;
      verbose && exit(verbose);
    } # for (chr in ...)
  } # for (aa in ...)
  verbose && exit(verbose);

  res;
}, protected=TRUE)  # calculateChromosomeStatistics()


##############################################################################
# HISTORY:
# 2007-09-16
# o Now process() of CopyNumberSegmentationModel reports timing information
#   for each chromosome fitted.
# o Made a better job cleaning out non-needed objects in process().
# 2007-09-15
# o Added getGenome() and setGenome().  Now, getGenomeData() searches for a
#   matching <genome>(,<tags)*,chromosomes.txt in 
#   annotationData/genomes/<genome>/. As a backup it also searches in the
#   corrsponding package directory(ies).
#   hgChromosomes.txt is no longer used.
# o BUG FIX: getGenomeData() was declared static, but it wasn't used that way.
# 2007-09-05
# o Was thinking to add a default asterisk tag to the output data set name.
#   However, although this works beautifully, the care has to be taken to
#   redesign ChromosomeExplorer, e.g. do you want to treat GLAD and CBS data
#   as totally different data sets?!? Possibly, because it is more consistent
#   with everything else, but for now we leave it as it since that works well.
# 2007-09-04
# o Now plot() is fully implemented CopyNumberSegmentationModel.  Subclasses
#   pretty much only have to implement pointsRawCNs() if wanted extra 
#   features, e.g. colors and outliers.
# 2007-08-20
# o It was actually not much that was hardwired to the GLAD model.
#   Note that most methods, including fit(), are generic enough to be 
#   defined in the superclass.  Methods that need to be implemented in
#   subclasses are: fitOne().
# o Created from GladModel.R.
##############################################################################
