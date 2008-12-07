setMethodS3("justSNPRMA", "character", function(...) {
  oligo::justSNPRMA(...);
})


setMethodS3("justSNPRMA", "AffymetrixCelSet", function(this, ..., normalizeToHapmap=TRUE, returnESet=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'normalizeToHapmap':
  normalizeToHapmap <- Arguments$getLogical(normalizeToHapmap);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Running SNPRMA on ", class(this)[1]);

  csR <- this;
  cdf <- getCdf(csR);
  chipType <- getChipType(cdf);
  hasCNs <- (regexpr("^GenomeWideSNP_(5|6)$", chipType) != -1);

  # Get the SNP only tag
  if (hasCNs) {
    snpOnlyTag <- "SNPs";
  } else {
    snpOnlyTag <- NULL;
  }
  snpOnlyTag <- "SNPs";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Rank-based quantile normalization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Rank-based quantile normalization");

  verbose && enter(verbose, "Setting up normalization model");
  if (normalizeToHapmap) {
    refTag <- "HapMapRef";
  } else {
    refTag <- NULL;
  }
  qn <- QuantileNormalization(csR, targetDistribution=NULL,
                                subsetToAvg=NULL, typesToUpdate="pm",
                                        tags=c("*", snpOnlyTag, refTag));

  if (!isDone(qn)) {
    if (normalizeToHapmap) {
      verbose && enter(verbose, "Loading HapMap reference target quantiles");
  
      pdPkgName <- oligo::cleanPlatformName(chipType);
      verbose && cat(verbose, "Platform Design (PD) package: ", pdPkgName);
  
      # Load target from PD package
      path <- system.file(package=pdPkgName);
      if (path == "") {
        throw("Cannot load HapMap reference target quantiles. Package not installed: ", pdPkgName);
      }

      path <- file.path(path, "extdata");
      path <- Arguments$getReadablePath(path);
  
      verbose && enter(verbose, "Loading binary file");
      filename <- sprintf("%sRef.rda", pdPkgName);
      pathname <- Arguments$getReadablePathname(filename, path=path);
      verbose && cat(verbose, "Pathname: ", pathname);
      target <- loadToEnv(pathname)$reference;
      verbose && str(verbose, target);
      verbose && exit(verbose);

      qn$.targetDistribution <- target;
      rm(target);
      verbose && exit(verbose);
    }

    verbose && enter(verbose, "Identifying cells for fitting normalization function");

    # justSNPRMA() operates only on SNP_A-* units (e.g. CN units ignored).
    # For this reason we here *estimate* the normalization function based
    # on these units only, but for convenience we will apply it to all
    # units (including CN units, if they exist).
    verbose && enter(verbose, "Identifying units");
    pattern <- "^SNP_A-";
    verbose && cat(verbose, "Pattern: ", pattern);
    units <- indexOf(cdf, pattern=pattern);
    verbose && cat(verbose, "Units:");
    verbose && str(verbose, units);
    verbose && exit(verbose);

    verbose && enter(verbose, "Identifying cell indices of these units");
    cells <- getCellIndices(cdf, units=units, unlist=TRUE, useNames=FALSE);
    verbose && cat(verbose, "Cells:");
    verbose && str(verbose, cells);
    rm(units);
    verbose && exit(verbose);

    qn$.subsetToAvg <- cells;
    rm(cells);
    verbose && exit(verbose);
  }

  verbose && print(verbose, qn);
  verbose && exit(verbose);

  verbose && enter(verbose, "Processing");
  csN <- process(qn, verbose=verbose);
  verbose && print(verbose, csN);
  verbose && exit(verbose);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Probe-level summarization
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Fitting probe-level (summarization) model");
  # We use the oligo estimator for fitting the log-additive model
  plm <- RmaSnpPlm(csN, mergeStrands=FALSE, flavor="oligo");
  verbose && print(verbose, plm);

  if (length(findUnitsTodo(plm)) > 0) {
    if (hasCNs) {
      verbose && enter(verbose, "Fitting CN probes");
      units <- fitCnProbes(plm, verbose=verbose);
      verbose && cat(verbose, "CN units fitted:");
      verbose && str(verbose, units);
      verbose && exit(verbose);
    }
  
    verbose && enter(verbose, "Fitting remaining units");
    units <- fit(plm, verbose=verbose);
    verbose && cat(verbose, "Units fitted:");
    verbose && str(verbose, units);
    verbose && exit(verbose);
  }

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting chip effect set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Extracting ChipEffectSet");
  ces <- getChipEffectSet(plm);
  verbose && print(verbose, ces);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extracting eSet
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (returnESet) {
    verbose && enter(verbose, "Extracting eSet");
    if (hasCNs) {
      eSet <- extractSnpCnvQSet(ces, verbose=log);
    } else {
      eSet <- extractSnpQSet(ces, verbose=log);
    }
    verbose && print(verbose, eSet);
    verbose && exit(verbose);

    res <- eSet;
  } else {
    res <- ces;
  }

  # Return result
  res;
})


############################################################################
# HISTORY:
# 2008-12-05
# o Created.
############################################################################
