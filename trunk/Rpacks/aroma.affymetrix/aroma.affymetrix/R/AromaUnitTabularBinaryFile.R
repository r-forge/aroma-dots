###########################################################################/**
# @RdocClass AromaUnitTabularBinaryFile
#
# @title "The AromaUnitTabularBinaryFile class"
#
# \description{
#  @classhierarchy
#
#  A AromaUnitTabularBinaryFile is an @see "AromaTabularBinaryFile" with
#  the constraint that the rows map one-to-one to, and in the same order as,
#  the units in a CDF.  The (full) chip type of the CDF is given by the
#  mandatory file footer \code{chipType}.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaTabularBinaryFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \seealso{
#   @see "AromaTabularBinaryFile".
# }
#*/########################################################################### 
setConstructorS3("AromaUnitTabularBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), "AromaUnitTabularBinaryFile",
    "cached:.cdf" = NULL
  );
})


setMethodS3("as.character", "AromaUnitTabularBinaryFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Chip type: %s", getChipType(this)));
  n <- length(s);
  s <- s[c(1:(n-2), n, n-1)];

  class(s) <- class;
  s;
})


setMethodS3("clearCache", "AromaUnitTabularBinaryFile", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".cdf")) {
    this[[ff]] <- NULL;
  }

  # Then for this object
  NextMethod(generic="clearCache", object=this, ...); 
}, private=TRUE)


setMethodS3("getFilenameExtension", "AromaUnitTabularBinaryFile", static=TRUE, abstract=TRUE);

setMethodS3("nbrOfUnits", "AromaUnitTabularBinaryFile", function(this, ...) {
  nbrOfRows(this, ...);
})


setMethodS3("getChipType", "AromaUnitTabularBinaryFile", function(this, fullname=TRUE, .old=FALSE, ...) {
  footer <- readFooter(this);
  chipType <- footer$chipType;

  if (is.null(chipType) && !.old) {
    throw("File format error: This ", class(this)[1], " file does not contain information on chip type in the file footer.  This is because the file is of an older file format an is no longer supported.  Please update to a more recent version: ", getPathname(this));
  }

  if (.old) {
    # Keep backward compatible for a while. /HB 2008-01-19
    if (fullname) {
      # Currently these files do not contain fullname chiptype information.
      # Instead we have to search for a match CDF and use its chip type.
      # Note, the CDF returned will depend of which CDF exists in the search
      # path, if at all.  AD HOC! /HB 2007-12-10
      cdf <- getCdf(this, ...);
      chipType <- getChipType(cdf, fullname=fullname);
    } else {
      chipType <- getName(this, ...);
    }
    chipType <- trim(chipType);
    if (nchar(chipType) == 0)
      throw("File format error: The inferred chip type is empty.");
  } else {
    if (!fullname) {
      chipType <- gsub(",.*", "", chipType);
    }
    chipType <- trim(chipType);
    if (nchar(chipType) == 0)
      throw("File format error: The chip type according to the file footer is empty.");
  }

  chipType;
})



setMethodS3("getCdf", "AromaUnitTabularBinaryFile", function(this, ..., force=FALSE, .old=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  cdf <- this$.cdf;
  if (force || is.null(cdf)) {
    if (.old) {
      # Generate all possible fullname 'chipTypes' and search for the existance
      # of a CDF with the longest name *and* that have the same number of units.
      # This is AD HOC! We should really store the full chiptype of the
      # corresponding CDF in the internal file header. /HB 2007-12-10
  
      verbose && enter(verbose, "Searching for a match CDF");
  
      verbose && cat(verbose, "Filename: ", getFilename(this));
      name <- getName(this, ...);
      tags <- getTags(this, collapse=NULL, ...); 
  
      validator <- function(cdf, ...) {
        (nbrOfUnits(cdf) == nbrOfUnits(this));
      }
      pathname <- findByCdf2(chipType=name, tags=tags, validator=validator, 
                                                    verbose=less(verbose, 1));
      if (is.null(pathname)) {
        throw("Failed to locate a CDF for ", class(this)[1], 
              " that have ", nbrOfUnits, " units: ", getFullName(this));
      }
  
      cdf <- AffymetrixCdfFile$fromFile(pathname);
    
      verbose && exit(verbose);
    } else {
      cdf <- AffymetrixCdfFile$fromChipType(getChipType(this));
    }

    this$.cdf <- cdf;
  }

  cdf;
})



setMethodS3("fromChipType", "AromaUnitTabularBinaryFile", function(this, ...) {
  byChipType(this, ...);
})

setMethodS3("byChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, validate=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'chipType':
  chipType <- Arguments$getCharacter(chipType);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 

  pathnames <- findByChipType(static, chipType=chipType, tags=tags, firstOnly=FALSE, ...);
  if (is.null(pathnames)) {
    throw("Could not locate a file for this chip type: ", 
                                   paste(c(chipType, tags), collapse=","));
  }

  verbose && cat(verbose, "Number of tabular binary files located: ", length(pathnames));
  verbose && print(verbose, pathnames);

  # Validate?
  if (validate) {
    verbose && cat(verbose, "Validation against CDF requested");
    verbose && enter(verbose, "Locating all matching CDFs");

    # Locate corresponding CDFs
    fullname <- paste(c(chipType, tags), collapse=",");
    parts <- unlist(strsplit(fullname, split=",", fixed=TRUE));
    chipType <- parts[1];
    tags <- parts[-1];
    cdfPathnames <- findByCdf2(chipType, tags=tags, firstOnly=FALSE, 
                                                 verbose=less(verbose, 1));
    queryStr <- paste(c(chipType, tags), collapse=", ");
    if (is.null(cdfPathnames)) {
      throw("Cannot validate against CDF.  No CDF located: ", queryStr);
    }

    cdfs <- lapply(cdfPathnames, FUN=function(pathname) {
      AffymetrixCdfFile$fromFile(pathname, verbose=less(verbose, 5));
    });

    verbose && cat(verbose, "Number of CDFs found: ", length(cdfs));

    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Scanning for a valid file");

  for (kk in seq(along=pathnames)) {
    pathname <- pathnames[kk];
    verbose && enter(verbose, "File #", kk, " (", pathname, ")");

    # Create object
    res <- newInstance(static, pathname);

    # Validate?
    if (validate) {
      for (ll in seq(along=cdfs)) {
        cdf <- cdfs[[ll]];
        if (nbrOfUnits(res) == nbrOfUnits(cdf)) {
          verbose && cat(verbose, "Found a tabular binary file that is compatible with matching CDF:");
          verbose && cat(verbose, "Tabular binary file: ", getPathname(res));
          verbose && cat(verbose, "CDF: ", getPathname(cdf));
          break;
        }
        res <- NULL;
      }
    }

    if (!is.null(res)) {
      verbose && cat(verbose);
    }

    verbose && exit(verbose);
  }

  if (is.null(res)) {
    throw("Failed to located a (valid) tabular binary file: ", queryStr);
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)



setMethodS3("findByChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  # Strip 'monocell' parts
  parts <- parts[parts != "monocell"];
  chipType <- parts[1];
  tags <- parts[-1];
  fullname <- paste(c(chipType, tags), collapse=",");

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
#  args$firstOnly <- FALSE;
#  str(args);
  pathname <- do.call("findAnnotationDataByChipType", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    pattern <- sprintf("^%s.*[.]%s[.]lnk$", chipType, ext);
    args$pattern <- pattern;
    pathname <- do.call("findAnnotationDataByChipType", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- filePath(pathname, expandLinks="any");
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE)


setMethodS3("indexOfUnits", "AromaUnitTabularBinaryFile", function(this, names, ...) {
  # Look up unit names from CDF
  cdf <- getCdf(this);
  idxs <- match(names, getUnitNames(cdf));
  idxs;
}, protected=TRUE)




###########################################################################/**
# @RdocMethod allocateFromCdf
#
# @title "Creates an AromaUnitTabularBinaryFile mapping to a given CDF"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cdf}{The @see "AffymetrixCdfFile" used as a template and from
#      which the (full) chip type is taken.}
#   \item{path}{The path where to store the new file.}
#   \item{tags}{A @character @vector of optional tags appended to 
#      the filename.}
#   \item{footer}{A nested named @list structure of additional attributes
#      that are saved in the file footer after the mandatory ones.}
#   \item{...}{Additional arguments passed to \code{allocate()} of 
#      @see "AromaTabularBinaryFile".}
# }
#
# \value{
#  Returns a @see "AromaUnitTabularBinaryFile" object.
# }
#
# @author
#
# \seealso{
#   To update to file footer afterwards, see @seemethod "writeFooter".
#   @seeclass
# }
#
# @keyword IO
#*/########################################################################### 
setMethodS3("allocateFromCdf", "AromaUnitTabularBinaryFile", function(static, cdf, path=getPath(cdf), tags=NULL, footer=list(), ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(cdf);
  # Exclude 'monocell' tags
  chipType <- gsub(",monocell", "", chipType);

  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create tabular binary file
  res <- allocate(static, filename=filename, path=path, 
                                           nbrOfRows=nbrOfUnits(cdf), ...);


  # Write attributes to footer
  attrs <- list(chipType=chipType);
  footer <- c(attrs, footer);
  writeFooter(res, footer);

  res;
}, static=TRUE)




setMethodS3("importFrom", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  if (inherits(src, "AffymetrixNetAffxCsvFile")) {
    importFromAffymetrixNetAffxCsvFile(this, src, ...);
  } else if (inherits(src, "DChipGenomeInformation")) {
    importFromDChipGenomeInformation(this, src, ...);
  } else if (inherits(src, "GenomeInformation")) {
    importFromGenomeInformation(this, src, ...);
  } else if (inherits(src, "AffymetrixTabularFile")) {
    importFromAffymetrixTabularFile(this, src, ...);
  } else if (inherits(src, "GenericTabularFile")) {
    importFromGenericTabularFile(this, src, ...);
  } else {
    throw("Do not know how to import from an src of class ", class(src)[1]);
  }
})


setMethodS3("importFromGenericTabularFile", "AromaUnitTabularBinaryFile", abstract=TRUE);


setMethodS3("importFromAffymetrixTabularFile", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  if (!inherits(src, "AffymetrixTabularFile")) {
    throw("Argument 'src' is not a AffymetrixTabularFile file: ", class(src)[1]);
  }

  importFromGenomeInformation(this, src, ...);
});


setMethodS3("importFromAffymetrixNetAffxCsvFile", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);

setMethodS3("importFromDChipGenomeInformation", "AromaUnitTabularBinaryFile", function(this, src, ...) {
  # Argument 'src':
  if (!inherits(src, "DChipGenomeInformation")) {
    throw("Argument 'src' is not a DChipGenomeInformation file: ", class(src)[1]);
  }

  importFromGenomeInformation(this, src, ...);
})


setMethodS3("importFromGenomeInformation", "AromaUnitTabularBinaryFile", abstract=TRUE, protected=TRUE);


############################################################################
# HISTORY:
# 2008-02-13
# o Added and updated Rdoc comments.
# 2008-01-19
# o Now AromaUnitTabularBinaryFile gets the chip type from the file footer.
# o ROBUSTNESS: Now fromChipType() of AromaUnitTabularBinaryFile validates
#   that the number of units in the located file match the number of units
#   in the CDF located using the same search parameters.
# 2007-12-10
# o Currently a AromaUnitTabularBinaryFile (e.g. AromaUgpFile) does not
#   contain information about the "fullname" chip type, but only the basic
#   chip-type name, e.g. we cannot infer the full chip-type name from 
#   'GenomeWideSNP_5,Full,r2.ugp', but only 'GenomeWideSNP_5'. The fullname
#   should be the same as the full chip-type name of the CDF used to define
#   the the unit map, e.g. 'GenomeWideSNP_5,Full.CDF'.
#   We should add a header (or footer) field in the file format that 
#   indicates the full chip type.  
#   However, until that is done, the best we can do is to turn to the ad
#   hoc solution of scanning for the CDF with the longest matching fullname,
#   if both 'GenomeWideSNP_5,Full.CDF' and 'GenomeWideSNP_5.CDF' exists,
#   the we match the former to 'GenomeWideSNP_5,Full,r2.ugp'.  The fullname
#   chip type of the UGP is then full chip-type name of the CDF.  NOTE,
#   there is major drawback with this.  If the user deletes the "full" CDF,
#   the above approach would all of a sudden return a different full name!
# o Added clearCache().
# 2007-09-14
# o Renames createFromCdf() to allocateFromCdf().
# 2007-09-13
# o Created from AromaUflFile.R.
############################################################################
