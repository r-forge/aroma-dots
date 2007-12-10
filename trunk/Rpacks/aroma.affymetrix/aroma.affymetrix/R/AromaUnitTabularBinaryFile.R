setConstructorS3("AromaUnitTabularBinaryFile", function(...) {
  extend(AromaTabularBinaryFile(...), "AromaUnitTabularBinaryFile",
    "cached:.cdf" = NULL
  );
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


setMethodS3("getChipType", "AromaUnitTabularBinaryFile", function(this, fullname=TRUE, ...) {
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

  chipType;
});


setMethodS3("getCdf", "AromaUnitTabularBinaryFile", function(this, ..., force=FALSE, verbose=FALSE) {
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
    # Generate all possible fullname 'chipTypes' and search for the existance
    # of a CDF with the longest name *and* that have the same number of units.
    # This is AD HOC! We should really store the full chiptype of the
    # corresponding CDF in the internal file header. /HB 2007-12-10
    verbose && enter(verbose, "Searching for a match CDF");
    verbose && cat(verbose, "Filename: ", getFilename(this));
    name <- getName(this, ...);
    tags <- getTags(this, collapse=NULL, ...);
    cdf <- NULL;
    # Number of units the matching CDF must have.
    nbrOfUnits <- nbrOfUnits(this);
    for (kk in rev(c(0,seq(along=tags)))) {
      cdfTags <- tags[seq(length=kk)];
      fullname <- paste(c(name, cdfTags), collapse=",");
      verbose && printf(verbose, "Trying '%s'...", fullname);
      tryCatch({
        cdf <- AffymetrixCdfFile$fromChipType(name, tags=cdfTags);
        # Does it have the correct number of units?
        if (nbrOfUnits(cdf) != nbrOfUnits)
          cdf <- NULL;
      }, error = function(ex) {})
  
      # Found a CDF?
      if (!is.null(cdf)) {
        verbose && writeRaw(verbose, "found.\n");
        break;
      }
      verbose && writeRaw(verbose, "no match.\n");
    }
    verbose && exit(verbose);
  
    if (is.null(cdf)) {
      throw("Failed to locate a CDF for ", class(this)[1], 
            " that have ", nbrOfUnits, " units: ", getFullName(this));
    }
  
    this$.cdf <- cdf;
  }

  cdf;
})


setMethodS3("fromChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  pathname <- findByChipType(static, chipType=chipType, tags=tags, ...);
  if (is.null(pathname)) {
    throw("Could not locate file for this chip type: ", 
                                   paste(c(chipType, tags), collapse=","));
  }

  # Create object
  newInstance(static, pathname);
}, static=TRUE)



setMethodS3("findByChipType", "AromaUnitTabularBinaryFile", function(static, chipType, tags=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/chipTypes/<chipType>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get fullname, name, and tags
  fullname <- paste(c(chipType, tags), collapse=",");
  parts <- unlist(strsplit(fullname, split=","));
  chipType <- parts[1];
  tags <- parts[-1];

  ext <- getFilenameExtension(static);
  ext <- paste(c(tolower(ext), toupper(ext)), collapse="|");
  ext <- sprintf("(%s)", ext);

  pattern <- sprintf("^%s.*[.]%s$", fullname, ext);
  args <- list(chipType=chipType, ...);
  args$pattern <- pattern;  # Override argument 'pattern'?
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



setMethodS3("allocateFromCdf", "AromaUnitTabularBinaryFile", function(static, cdf, path=getPath(cdf), tags=NULL, ...) {
  # Argument 'cdf':
  if (!inherits(cdf, "AffymetrixCdfFile")) {
    throw("Argument 'cdf' is not an AffymetrixCdfFile: ", class(cdf)[1]);
  }

  # Generate filename: <chipType>(,tags)*.<ext>
  chipType <- getChipType(cdf);
  fullname <- paste(c(chipType, tags), collapse=",");
  ext <- getFilenameExtension(static);
  filename <- sprintf("%s.%s", fullname, ext);

  # Create tabular binary file
  allocate(static, filename=filename, path=path, nbrOfRows=nbrOfUnits(cdf), ...);
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
