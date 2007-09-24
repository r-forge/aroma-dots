setMethodS3("getXTheta", "CopyNumberSegmentationModel", function(this, chromosome, reorder=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Get (position, chipType, unit) map
  pcu <- getPositionChipTypeUnit(this, chromosome=chromosome, 
                                                  verbose=less(verbose, 10));

  # Get list of chip-effect sets
  cesList <- getListOfChipEffects(this);

  # Allocate return structure
  x <- vector("double", nrow(pcu));
  theta <- matrix(NA, nrow=length(xs), ncol=nbrOfArrays(this));

  for (kk in seq(along=cesList)) {
    verbose && enter(verbose, "Chip type #", kk, " of ", length(cesList));
    ces <- cesList[[kk]];
    idxs <- which(as.integer(pcu[,"chipType"]) == kk);
    x[idxs] <- pcu[idxs,"position"];
    units <- pcu[idxs,"unit"];
    verbose && cat(verbose, "Units: ");
    verbose && str(verbose, units);

    verbose && enter(verbose, "Reading data across chip-effect files");
    theta[idxs,] <- extractMatrix(ces, units=units, verbose=less(verbose, 10));
    verbose && exit(verbose);

    verbose && exit(verbose);
  }
  verbose && exit(verbose);

  list(x=x, theta=theta);
}, protected=TRUE);  # getXTheta()


##############################################################################
# HISTORY:
# 2007-09-24
# o Added getXTheta().
##############################################################################
