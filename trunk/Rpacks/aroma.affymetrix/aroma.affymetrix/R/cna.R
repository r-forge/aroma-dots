###########################################################################/**
# @set "class=list"
# @RdocMethod cna
#
# @title "Copy-number analysis"
#
# \description{
#  This high-level function takes a list of probe-level data set and 
#  performs a complete copy-number analysis.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSets}{A @list of @see "AffymetrixCelSet".}
#   \item{pres}{A @character @vector describing what methods should be used
#     to preprocess the probe-level data.
#     If @NULL, this step is skipped.}
#   \item{plm}{A @character string describing what probe-level model (PLM)
#     to fit to the preprocessed data.}
#   \item{posts}{A @character @vector describing what methods should be used
#     to postprocess the chip-effect estimates of the PLM.
#     If @NULL, this step is skipped.}
#   \item{cnm}{A @character string describing what copy-number model (CNM)
#     to fit to identify copy-number regions.}
#   \item{...}{Not used.}
#   \item{ram}{A @double in (0,Inf) specifying if more or less (RAM) 
#     memory than default should be used.  The smaller the value is, 
#     the more memory conservative the analysis will be.}
#   \item{verbose}{Specifies how much details should be outputted while
#     analysing data.  The lower threshold, the more details.}
# }
#
# @author
#*/###########################################################################
setMethodS3("cna", "list", function(dataSets, pres="QN", plm=c("RMA", "MBEI"), posts="FLN", cnm="glad", ..., ram=1, verbose=-2) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSets':
  for (dataSet in dataSets) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSets' contains a non-AffymetrixCelSet: ", 
                                                           class(dataSet)[1]);
    }
  }

  # Argument 'pres':
  if (!is.null(pres)) {
    pres <- Arguments$getCharacters(pres);
  }

  # Argument 'plm':
  plm <- match.arg(plm);

  # Argument 'post':
  if (!is.null(posts)) {
    posts <- Arguments$getCharacters(posts);
  }

  # Argument 'cnm':
  cnm <- Arguments$getCharacter(cnm);

  # Argument 'ram':
  ram <- Arguments$getDouble(ram, range=c(0.001, Inf));

  # Argument 'verbose':
  if (inherits(verbose, "Verbose")) {
  } else {
    verbose <- Arguments$getVerbose(verbose);
    setTimestampOn(verbose);
  }
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  gc();


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Pre-processing (of probe-level data)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Preprocessing (probe-level data)");
  listOfProbeData <- dataSets;
  for (kk in seq(along=dataSets)) {
    probeData <- listOfProbeData[[kk]];
    for (pre in pres) {
      if (pre == "ACT") {
        model <- AllelicCrosstalkNormalization(probeData);
      } else if (pre == "QN") {
        model <- QuantileNormalization(probeData);
      }
      print(model);
      probeData <- process(model, verbose=verbose);
      rm(model);
      gc();
    }
    listOfProbeData[[kk]] <- probeData;
    rm(probeData);
    gc();
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Probe-level modelling
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Probe-level modelling");
  listOfCes <- list();
  for (kk in seq(along=dataSets)) {
    probeData <- listOfProbeData[[kk]];
    if (plm == "RMA") {
      Plm <- RmaCnPlm;
    } else if (plm == "MBEI") {
      Plm <- MbeiCnPlm;
    }
    model <- Plm(probeData, combineAlleles=TRUE, mergeStrands=TRUE);
    rm(probeData);
    listOfProbeData[[kk]] <- NULL;
    gc();
    print(model);
    fit(model, moreUnits=ram, verbose=verbose);
    gc();
    listOfCes[[kk]] <- getChipEffects(model);
    rm(model);
    gc();
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Post-processing (of chip-effect estimates)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Postprocessing (chip-effect estimates)");
  for (kk in seq(along=dataSets)) {
    ces <- listOfCes[[kk]];

    for (post in posts) {
      if (post == "FLN") {
        model <- FragmentLengthNormalization(ces);
      }
      print(model);
      ces <- process(model, verbose=verbose);
      rm(model);
      gc();
    }

    # Ad hoc for now
    ces$combineAlleles <- TRUE;
    ces$mergeStrands <- TRUE;

    listOfCes[[kk]] <- ces;
    rm(ces);
    gc();
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy-number modelling
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Identifying copy-number regions");
  if (cnm == "glad") {
    model <- MultiGladModel(listOfCes);
  }
  rm(listOfCes);
  gc();
  print(model);
  fit(model, verbose=verbose);
  gc();
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # ChromosomeExplorer
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Generating report");
  ce <- ChromosomeExplorer(model);
  rm(model);
  gc();
  print(ce);
  process(ce, verbose=verbose);
  verbose && exit(verbose);

  invisible(dataSets);
}) # cna()


setMethodS3("cna", "AffymetrixCelSet", function(dataSet, ...) {
  cna(list(dataSet), ...);
})


############################################################################
# HISTORY:
# 2007-01-04
# o Created.
############################################################################
