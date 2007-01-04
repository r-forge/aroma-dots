library(aroma.affymetrix);
source("init.R");
source("ui.R");

###########################################################################/**
# @RdocFunction cna
#
# @title "Copy-number analysis"
#
# \description{
#  This function takes a raw data set and performs a complete copy-number
#  analysis.
# }
# 
# @synopsis
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#   \item{...}{Not used.}
# }
#
# @author
# 
# \seealso{
#   An object of this class is typically part of a @see "CnChipEffectSet".
# }
#*/###########################################################################
cna <- function(dataSets=NULL, pres="qn", plm=c("rma","mbei"), posts="fln", cnm="glad", ..., ram=1, verbose=-2) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSets':
  if (is.null(dataSets)) {
    dataSets <- selectDataSets(c("raw", "probeData"));
  } else if (is.character(dataSets)) {
    # A path to a data set
    path <- dataSets;
    if (isDirectory(path)) {
      dataSets <- list(AffymetrixCelSet$fromFiles(path));
    }
  } else if (is.list(dataSets)) {
    dataSets <- list(dataSets);
  }

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
  verbose <- Arguments$getVerbose(verbose);
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
      if (pre == "act") {
        model <- AllelicCrosstalkNormalization(probeData);
      } else if (pre == "qn") {
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
    if (plm == "rma") {
      model <- RmaCnPlm(probeData, combineAlleles=TRUE, mergeStrands=TRUE);
    } else if (plm == "mbei") {
      model <- MbeiCnPlm(probeData, combineAlleles=TRUE, mergeStrands=TRUE);
    }
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
      if (post == "fln") {
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
} # cna()


############################################################################
# HISTORY:
# 2007-01-04
# o Created.
############################################################################
