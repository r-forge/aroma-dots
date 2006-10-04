###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractSnpFeatureSet
#
# @title "Extracts an in-memory SnpFeatureSet object from the data set"
#
# \description{
#  @get "title".
#  Note that any modifications done to the extract object will \emph{not}
#  be reflected in the original data set.
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

  verbose && enter(verbose, "Extracting SnpFeatureSet from data set");
  cdf <- getCdf(this);
  chipType <- getChipType(cdf);

  verbose && printf(verbose, "Chip type: %s\n", chipType);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert that the corresponding PD package is installed, e.g.
  # pdmapping50kxba240. 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Source: http://www.biostat.jhsph.edu/~bcarvalh/old/research.html
  safeChipType <- gsub("cdf$", "", cleancdfname(chipType));
  pdName <- sprintf("pd%s", safeChipType);
  verbose && enter(verbose, "Trying required probe-definition package: ", pdName);
  pkg <- packageDescription(pdName);
  if (identical(pkg, NA))
    throw("Probe-design package required by 'oligo' not installed: ", pdName);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Load PD environment
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # For a 100K chip this loads up ~700Mb of the memory.
  require(pdName, character.only=TRUE);
  gc();

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read data
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading SnpFeatureSet");
  # Work around for bug in 'oligo' not accepting filenames with a path.
  opwd <- setwd(getPath(this));
  on.exit(setwd(opwd), add=TRUE);
  filenames <- basename(getPathnames(this));
  res <- oligo::read.celfiles(filenames, ..., verbose=as.logical(less(verbose)));
  verbose && exit(verbose);

  verbose && exit(verbose);

  res;
})

############################################################################
# HISTORY:
# 2006-10-02
# o Wow, after quite a few hours troubleshooting on different systems it
#   turns out that read.celfiles() of oligo (v 0.99.20) does not work with
#   filenames with a path.
# o Created.
############################################################################
