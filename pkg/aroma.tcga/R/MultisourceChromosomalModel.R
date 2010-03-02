###########################################################################/**
# @RdocClass MultisourceChromosomalModel
#
# @title "The MultisourceChromosomalModel class"
#
# \description{
#  @classhierarchy
#
#  This \emph{abstract} class represents a chromosomal model based on
#  data from multiple sources (platforms, labs, methods, ...).
# }
# 
# @synopsis
#
# \arguments{
#   \item{dsList}{A @list of @see "AromaUnitSignalBinarySet".}
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
setConstructorS3("MultisourceChromosomalModel", function(dsList=NULL, tags="*", genome="Human", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dsList':
  if (!is.null(dsList)) {
    if (!is.list(dsList)) {
      throw("Argument 'dsList' must be a list: ", class(dsList)[1]);
    }

    className <- "AromaUnitSignalBinarySet";
    for (kk in seq(along=dsList)) {
      ds <- dsList[[kk]];
      if (!inherits(ds, className)) {
        throw(sprintf("Element #%d of argument 'dsList' is not a %s: %s", 
                                       kk, className, class(ds)[1]));
      }
    }
  }

  # Argument 'tags':
  tags <- Arguments$getTags(tags);
 

  this <- extend(Object(), "MultisourceChromosomalModel",
    .alias = NULL,
    .dsList = dsList,
    .chromosomes = NULL,
    .tags = tags,
    .genome = genome
  );

  # Validate?
  if (!is.null(this$.dsList)) {
    # Validate genome
    pathname <- getGenomeFile(this);
  }

  this;
}, abstract=TRUE)


setMethodS3("as.character", "MultisourceChromosomalModel", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, paste("Name:", getName(this)));
  s <- c(s, paste("Tags:", getTags(this, collapse=",")));
  s <- c(s, paste("Chip type (merged):", getChipType(this, merge=TRUE)));
  s <- c(s, sprintf("Path: %s", getPath(this)));

  s <- c(s, sprintf("Number of data sets (sources): %s", nbrOfSources(this)));
  s <- c(s, "List of data sets:");
  dsList <- getDataSetList(this);
  s <- c(s, sapply(dsList, as.character));

  s <- c(s, sprintf("RAM: %.2fMB", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";

  s;
}, protected=TRUE)


setMethodS3("clearCache", "MultisourceChromosomalModel", function(this, ...) {
  # Clear all cached values.
  # /AD HOC. clearCache() in Object should be enough! /HB 2007-01-16
  for (ff in c()) {
    this[[ff]] <- NULL;
  }

  if (!is.null(this$.dsList)) {
    clearCache(this$.dsList);
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...);
})


setMethodS3("getDataSetList", "MultisourceChromosomalModel", function(this, ...) {
  this$.dsList;
})


setMethodS3("nbrOfSources", "MultisourceChromosomalModel", function(this, ...) {
  dsList <- getDataSetList(this);
  length(dsList);
})


setMethodS3("getRootPath", "MultisourceChromosomalModel", function(this, ...) {
  tag <- getAsteriskTags(this)[1];
  sprintf("%sData", tolower(tag));
})


setMethodS3("getParentPath", "MultisourceChromosomalModel", function(this, ...) {
  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # The full path
  path <- filePath(rootPath, fullname, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path)) {
      throw("Failed to create directory: ", path);
    }
  }

  path;
})

setMethodS3("getPath", "MultisourceChromosomalModel", function(this, ...) {
  path <- getParentPath(this, ...);

  # Chip type
  chipType <- getChipType(this);

  # The full path
  path <- filePath(path, chipType, expandLinks="any");

  # Create path?
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path)) {
      throw("Failed to create output directory: ", path);
    }
  }

  path;
})

setMethodS3("getReportPath", "MultisourceChromosomalModel", function(this, ...) {
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



setMethodS3("getChipTypes", "MultisourceChromosomalModel", function(this, ...) {
  dsList <- getDataSetList(this);
  chipTypes <- sapply(dsList, FUN=getChipType, ...);

  chipTypes;
})

setMethodS3("getChipType", "MultisourceChromosomalModel", function(this, merge=TRUE, collapse="+", ...) {
  chipTypes <- getChipTypes(this, ...);

  # Merge to a single string?
  if (merge) {
    chipTypes <- mergeByCommonTails(chipTypes, collapse=collapse);
  }
 
  chipTypes;
})



###########################################################################/**
# @RdocMethod getNames
#
# @title "Gets the names of the arrays"
#
# \description{
#  @get "title" available to the model.
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
setMethodS3("getNames", "MultisourceChromosomalModel", function(this, ...) {
  dsList <- getDataSetList(this);
  names <- lapply(dsList, FUN=getNames);
  names <- unlist(names, use.names=FALSE);
  names <- unique(names);
  names <- sort(names);
  names;
})

setMethodS3("getFullNames", "MultisourceChromosomalModel", function(this, ...) {
  getNames(this, ...);
})

setMethodS3("indexOf", "MultisourceChromosomalModel", function(this, patterns=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getFullNames <- function(fullnames=NULL, ...) {
    if (!is.null(fullnames))
      return(fullnames);
    getFullNames(this);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  names <- getNames(this);

  # Return all indices
  if (is.null(patterns)) {
    res <- seq(along=names);
    names(res) <- names;
    return(res);
  } else if (is.numeric(patterns)) {
    n <- length(names);
    res <- Arguments$getIndices(patterns, max=n);
    names(res) <- names[res];
    return(res);
  }

  fullnames <- NULL;

  naValue <- as.integer(NA);

  patterns0 <- patterns;
  res <- lapply(patterns, FUN=function(pattern) {
    pattern <- sprintf("^%s$", pattern);
    pattern <- gsub("\\^\\^", "^", pattern);
    pattern <- gsub("\\$\\$", "$", pattern);

    # Specifying tags?
    if (regexpr(",", pattern) != -1) {
      fullnames <- getFullNames(fullnames);
      idxs <- grep(pattern, fullnames);
    } else {
      idxs <- grep(pattern, names);
    }
    if (length(idxs) == 0)
      idxs <- naValue;

    # Note that 'idxs' may return more than one match
    idxs;
  });

  ns <- sapply(res, FUN=length);
  names <- NULL;
  for (kk in seq(along=ns)) {
    names <- c(names, rep(patterns0[kk], times=ns[kk]));
  }
  res <- unlist(res, use.names=FALSE);
  names(res) <- names;

  # Not allowing missing values?
  if (any(is.na(res))) {
    names <- names(res)[is.na(res)];
    throw("Some names where not match: ", paste(names, collapse=", "));
  }

  res;
})




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
setMethodS3("nbrOfArrays", "MultisourceChromosomalModel", function(this, ...) {
  length(getNames(this, ...));
})


setMethodS3("setName", "MultisourceChromosomalModel", function(this, name, ...) {
  # Argument 'name':
  name <- Arguments$getCharacter(name);
  this$.name <- name;
  invisible(this);
})


setMethodS3("getName", "MultisourceChromosomalModel", function(this, collapse="+", ...) {
  name <- this$.name;
  if (is.null(name)) {
    dsList <- getDataSetList(this);
    ds <- dsList[[1]];
    name <- getName(ds, ...);
  }
  name;
})

setMethodS3("getAsteriskTags", "MultisourceChromosomalModel", function(this, collapse=NULL, ...) {
  # Create a default asterisk tags for any class by extracting all
  # capital letters and pasting them together, e.g. AbcDefGhi => ADG.
  name <- class(this)[1];

  # Remove any 'Model' suffixes
  name <- gsub("Model$", "", name);

  name <- capitalize(name);

  # Vectorize
  name <- strsplit(name, split="")[[1]];

  # Identify upper case
  name <- name[(toupper(name) == name)];

  # Paste
  name <- paste(name, collapse="");

  name;
}, protected=TRUE)




setMethodS3("getTags", "MultisourceChromosomalModel", function(this, collapse=NULL, ...) {
  dsList <- getDataSetList(this);
  tags <- lapply(dsList, FUN=getTags, collapse=NULL);

  # Keep common tags
  tags <- getCommonListElements(tags);
  tags <- tags[[1]];

  # Add model tags
  tags <- c(tags, this$.tags);

  # Update default tags
  asteriskTags <- getAsteriskTags(this, collapse=",");
  if (length(asteriskTags) == 0)
    asteriskTags <- "";
  tags[tags == "*"] <- asteriskTags;

  tags <- Arguments$getTags(tags, collapse=NULL);

  # Get unique tags
  tags <- locallyUnique(tags);

  # Collapsed?
  tags <- Arguments$getTags(tags, collapse=collapse);

  tags;
})


setMethodS3("getFullName", "MultisourceChromosomalModel", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
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
setMethodS3("getChromosomes", "MultisourceChromosomalModel", function(this, ...) {
  dsList <- getDataSetList(this);
  ugpList <- lapply(dsList, FUN=getAromaUgpFile);
  chromosomes <- lapply(ugpList, getChromosomes);
  chromosomes <- unlist(chromosomes, use.names=TRUE);
  chromosomes <- sort(unique(chromosomes));
  chromosomes;
})


setMethodS3("getGenome", "MultisourceChromosomalModel", function(this, ...) {
  this$.genome;
})


setMethodS3("getGenomeFile", "MultisourceChromosomalModel", function(this, ..., verbose=FALSE) {
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


setMethodS3("setGenome", "MultisourceChromosomalModel", function(this, genome, tags=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'genome':
  genome <- Arguments$getCharacter(genome, length=c(1,1));

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



setMethodS3("getGenomeData", "MultisourceChromosomalModel", function(this, ..., verbose=FALSE) {
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


setMethodS3("fit", "MultisourceChromosomalModel", abstract=TRUE);


setMethodS3("getSetTag", "MultisourceChromosomalModel", function(this, ...) {
  tolower(getAsteriskTags(this)[1]);
}, private=TRUE)




setMethodS3("extractRawGenomicSignals", "MultisourceChromosomalModel", function(this, name, chromosome, region=NULL, ..., extractFcn=extractRawGenomicSignals, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'name':
  name <- Arguments$getCharacter(name);

  # Argument 'chromosome':
  allChromosomes <- getChromosomes(this);
  chromosome <- Arguments$getIndex(chromosome, range=range(allChromosomes));
  if (!is.element(chromosome, allChromosomes)) {
    throw("The value of argument 'chromosome' is unknown: ", chromosome);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Extract raw genomic signals");
  dsList <- getDataSetList(this);
  dsList <- lapply(dsList, FUN=extract, name);
  dfList <- lapply(dsList, FUN=getFile, 1);

  cnList <- NULL;
  for (kk in seq(along=dfList)) {
    df <- dfList[[kk]];
    if (!isFile(df))
      next;

    cn <- extractFcn(df, chromosome=chromosome, region=region);
    cn$id <- rep(kk, nbrOfLoci(cn));
    cnList <- c(cnList, list(cn));
  }

  res <- Reduce(append, cnList);

  verbose && print(verbose, res);
  verbose && exit(verbose);

  res;
}, protected=TRUE)



##############################################################################
# HISTORY:
# 2010-01-25
# o Created from ChromosomalModel.
##############################################################################
