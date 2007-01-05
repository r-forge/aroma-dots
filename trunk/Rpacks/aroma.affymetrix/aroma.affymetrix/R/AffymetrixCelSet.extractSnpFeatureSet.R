###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractSnpFeatureSet
#
# @title "Extracts an in-memory SnpFeatureSet object from the CEL set"
#
# \description{
#  @get "title".
#  Note that any modifications done to the extract object will \emph{not}
#  be reflected in the original CEL set.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Argument passed to \code{oligo::read.cellfiles()}.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "oligo::SnpFeatureSet-class" object.
# }
#
# \details{
#  Since the \pkg{oligo} package is making use of special CDF environment
#  packages, this method will warn if the needed package is missing and
#  explain that \pkg{oligo} will later try to download and install it
#  automatically.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#
# @keyword IO
# @keyword programming
#*/###########################################################################
setMethodS3("extractSnpFeatureSet", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  require(affy) || throw("Package not loaded: affy");
  require(oligo) || throw("Package not loaded: oligo");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Extracting SnpFeatureSet from CEL set");
  cdf <- getCdf(this);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that the corresponding PD package is installed, e.g.
  # pdmapping50kxba240. 
  # Source: http://www.biostat.jhsph.edu/~bcarvalh/old/research.html
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Loading platform-design package");
  pd <- PlatformDesign(cdf);
  pdName <- getPackageName(pd);
  verbose && cat(verbose, "Package: ", pdName);
  if (!isInstalled(pd))
    throw("Probe-design package required by 'oligo' not installed: ", pdName);
  # Load platform-design package (a 100K chip requires ~700MB RAM).
  load(pd);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading SnpFeatureSet");
  # Work around for bug in 'oligo' not accepting filenames with a path.
  # Note that this requires that all CEL files are in the same directory!
  opwd <- setwd(getPath(this));
  on.exit(setwd(opwd), add=TRUE);
  filenames <- basename(getPathnames(this));
  verbose && cat(verbose, "Current path: ", getwd());
  verbose && cat(verbose, "Filenames: ", paste(filenames, collapse=", "));
  res <- oligo::read.celfiles(filenames, ..., verbose=as.logical(less(verbose)));
  verbose && exit(verbose);

  # Updating the sample names
  sampleNames(res) <- getNames(this);

  verbose && exit(verbose);

  res;
})

############################################################################
# HISTORY:
# 2006-12-12
# o Now extractSnpFeatureSet() sets the sample names (instead of using the
#   filenames).
# 2006-12-11
# o extractSnpFeatureSet() is now using of the PlatformDesign class.
# 2006-10-02
# o Wow, after quite a few hours troubleshooting on different systems it
#   turns out that read.celfiles() of oligo (v 0.99.20) does not work with
#   filenames with a path.
# o Created.
############################################################################
