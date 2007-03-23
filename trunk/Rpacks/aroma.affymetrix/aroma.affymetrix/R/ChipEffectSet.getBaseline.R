###########################################################################/**
# @set "class=ChipEffectSet"
# @RdocMethod getBaseline
#
# @title "Gets the baseline signals across chromosomes"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{force}{If @TRUE, the CEL file that stores the is recreated.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getBaseline", "ChipEffectSet", function(this, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting CEL file to store baseline signals");
  path <- getPath(this);
  key <- list(dataset=getFullName(this), samples=getNames(this));
  id <- digest(key);
  filename <- sprintf(".baseline,%s.cel", id);
  pathname <- Arguments$getWritablePathname(filename, path=path);

  # Get a template CEL file
  df <- getFile(this, 1);

  if (!force && isFile(pathname)) {
    verbose && enter(verbose, "Retrieving existing CEL file");
    res <- fromFile(df, filename=pathname, path=NULL, 
                                                  verbose=less(verbose));
    verbose && exit(verbose);
  } else {
    verbose && enter(verbose, "Creating empty CEL file");
    # Create new empty file from this one
    res <- createFrom(df, filename=pathname, path=NULL, 
                                                  verbose=less(verbose));
    verbose && exit(verbose);
  }
  verbose && print(verbose, res);
  rm(df);
  verbose && exit(verbose);

  res;
})



setMethodS3("getBaseline", "SnpChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline", this, ...);
  res$mergeStrands <- getMergeStrands(this);
  res;
})

setMethodS3("getBaseline", "CnChipEffectSet", function(this, ...) {
  res <- NextMethod("getBaseline", this, ...);
  res$combineAlleles <- getCombineAlleles(this);
  res;
})




############################################################################
# HISTORY:
# 2007-03-22
# o TO DO: Estimate standard errors just like getAverage() does.
# o Added getBaseline().
# o First working version of calculateBaseline(). Method now creates a CEL
#   files to store the estimates.
# 2007-03-16
# o Created.  See ploidy4.ps paper.
############################################################################
