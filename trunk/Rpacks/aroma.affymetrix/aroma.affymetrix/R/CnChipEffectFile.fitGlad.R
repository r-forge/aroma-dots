###########################################################################/**
# @set "class=CnChipEffectFile"
# @RdocMethod fitGlad
#
# @title "Fits the GLAD model to copy-number estimates"
#
# \description{
#  @get "title" of a certain chromosome.
# }
#
# @synopsis
#
# \arguments{
#   \item{reference}{The @see "CnChipEffectFile" object used as the 
#     reference chip effects when calculating the raw (relative) copy-number
#     estimates..}
#   \item{chromosomes}{The chromosomes for which the model should be fitted.}
#   \item{units}{(Optional) The subset of units to be matched.}
#   \item{...}{Not used.}
#   \item{force}{If @TRUE, any in-memory cached results are ignored.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns the @see "GLAD::profileCGH" object returned by @see "GLAD::glad".
# }
#
# @author
#
# \seealso{
#   Internally @see "GLAD::glad" is used.
#   @seeclass
# }
#*/###########################################################################
setMethodS3("fitGlad", "CnChipEffectFile", function(this, reference, chromosomes=c(1:22,"X"), units=NULL, useStddvs=TRUE, ..., force=FALSE, verbose=FALSE) {
  throw("fitGlad() for CnChipEffectFile is deprecated since 2006-12-15.  Use the GladModel class instead.");

  require(GLAD) || throw("Package 'GLAD' not loaded.");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (!inherits(reference, "CnChipEffectFile")) {
    throw("Argument 'reference' is not a CnChipEffectFile: ", 
                                                        class(reference)[1]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  if (length(chromosomes) > 1) {
    res <- list();
    verbose && enter(verbose, "Fitting ", length(chromosomes), " chromosomes");
    verbose && cat(verbose, "Chromosomes: ", paste(chromosomes, collapse=","));
    for (cc in chromosomes) {
      ccStr <- as.character(cc);
      res[[ccStr]] <- fitGlad(this, reference=reference, chromosomes=cc, 
                       units=units, ..., force=force, verbose=less(verbose));
    }
    verbose && exit(verbose);
    return(res);
  }

  verbose && enter(verbose, "Fitting GLAD");

  chromosome <- chromosomes[1];

  cdf <- getCdf(this);
  chipType <- getChipType(cdf);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Extract arguments for glad().
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list(...);
  keep <- (names(args) %in% names(formals(glad.profileCGH)));
  gladArgs <- args[keep];

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached values
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  fullName <- getFullName(this);
  key <- list(method="fitGlad", 
              class=class(this)[1],
              dataSet=basename(dirname(getPath(this))),
              fullname=fullName,
              refset=basename(dirname(getPath(reference))),
              reference=getFullName(reference),
              chipType=chipType,
              chromosome=chromosome,
              units=units
             );
  dirs <- c("aroma.affymetrix", fullName);
  fit <- loadCache(key=key, dirs=dirs);
  if (!is.null(fit) && !force) {
    verbose && cat(verbose, "Cached on file.");
    verbose && exit(verbose);
    return(fit);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Fit GLAD
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get (X, M, A)
  verbose && enter(verbose, "Retrieving relative copy-number estimates");
  df <- getXAM(this, other=reference, chromosome=chromosome, units=units, 
                                                      verbose=less(verbose));

  verbose && enter(verbose, "Retrieving stddvs of chip effects");
  units <- as.integer(rownames(df));
  sdTheta <- getDataFlat(this, units=units, fields="sdTheta", verbose=less(verbose))[,"sdTheta"];
  df <- cbind(df, sdTheta);
  rm(sdTheta);
  verbose && exit(verbose);

  verbose && enter(verbose, "Re-order by physical position");
  df <- df[order(df[,"x"]),];
  verbose && exit(verbose);
  verbose && str(verbose, df);
  verbose && cat(verbose, sprintf("Extracted data for %d SNPs", nrow(df)));

  # Put the data in a format recognized by GLAD
  df <- data.frame(
    LogRatio=unname(df[,"M"]), 
    PosOrder=1:nrow(df), 
    Chromosome=rep(chromosome, nrow(df)), 
    PosBase=unname(df[,"x"]),
    # Add (chipType, units) identifiers to be able to backtrack SNP IDs etc.
    chipType=as.factor(chipType),
    units=units,
    sdTheta=unname(df[,"sdTheta"])
  );
  verbose && str(verbose, df);
  df <- as.profileCGH(df);
  verbose && str(verbose, df);
  verbose && exit(verbose);

  verbose && enter(verbose, "Calling glad()");
  args <- c(list(df), gladArgs, list(verbose=as.logical(verbose)));
  fit <- do.call("glad", args);
  verbose && exit(verbose);

  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  comment <- paste(unlist(key), collapse=";");
  saveCache(fit, key=key, comment=comment, dirs=dirs);

  fit;  
}, private=TRUE) # fitGlad()


############################################################################
# HISTORY:
# 2006-12-15
# o Made fitGlad() for CnChipEffectFile deprecated.  Use the GladModel class
#   instead.
# 2006-11-06
# o Now fitGlad() accepts a vector of chromosomes, because due to file 
#   caching it is probably faster to fit all chromosomes on one file than
#   across files.
# 2006-10-30
# o Now chip type and (CDF) unit indices are stored in the result too so
#   that the SNP IDs etc can be backtracked.
# o Added Rdoc comments.
# o BUG FIX: glad() requires data points to be order by physical position,
#   which was not (necessarily) the case.
# o Added more argument validation.
# 2006-10-17
# o Created.
############################################################################
