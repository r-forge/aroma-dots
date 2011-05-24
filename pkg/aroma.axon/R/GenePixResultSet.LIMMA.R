setMethodS3("readLimmaRGList", "GenePixResultSet", function(this, source="genepix.median", ..., translate=TRUE, fullnames=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'translate':
  translate <- Arguments$getLogical(translate);

  # Argument 'fullnames':
  fullnames <- Arguments$getLogical(fullnames);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  
  verbose && enter(verbose, sprintf("Reading %s as a limma RGList object", class(this)[1]));
  verbose && print(verbose, this);

  names <- NULL;
  if (translate) {
    if (fullnames) {
      names <- getFullNames(this);
    } else {
      names <- getNames(this);
      if (any(duplicated(names))) {
        throw("Cannot use 'fullnames=FALSE' when there exist duplicated names: ", hpaste(names[duplicated(names)]));
#        names <- getFullNames(this);
      }
    }
  }

  # Extract signals from all arrays
  pathnames <- getPathnames(this);

  verbose && enter(verbose, "Calling limma::read.maimages()");
#  verbose && cat(verbose, "Filenames: ", hpaste(basename(pathnames)));
  verbose && cat(verbose, "Pathnames:");
  pathnames <- unname(pathnames);
  verbose && str(verbose, pathnames);

  verbose && cat(verbose, "Names: ", hpaste(names));

  verboseT <- as.logical(verbose);
  data <- limma::read.maimages(files=pathnames, source=source, names=names, verbose=verboseT, ...);
  verbose && exit(verbose);

  verbose && enter(verbose, "Trimming read data");
  genes <- data$genes;
  genes$ID <- trim(genes$ID);
  genes$ID <- gsub("N/A", NA, genes$ID, fixed=TRUE);
  genes$Name <- trim(genes$Name);
  data$genes <- genes;
  rm(genes);
  verbose && exit(verbose);

  verbose && printf(verbose, "Object size: %d bytes\n", objectSize(data));

  verbose && exit(verbose);

  data;
}, protected=TRUE) # readLimmaRGList()



############################################################################
# HISTORY:
# 2011-05-23
# o Now readLimmaRGList() trims the gene annotations.
# o Renamed readLimmaData() to readLimmaRGList().
# o Added readLimmaData() from old GenePixDataSet.R.
# o Created.
############################################################################
