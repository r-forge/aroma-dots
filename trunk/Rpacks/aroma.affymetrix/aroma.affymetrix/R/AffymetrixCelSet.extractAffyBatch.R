###########################################################################/**
# @RdocMethod extractAffyBatch
#
# @title "Extracts an in-memory AffyBatch object from the data set"
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
#   \item{...}{Argument passed to \code{ReadAffy()} 
#     (@see "affy::read.affybatch").}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an @see "affy::AffyBatch-class" object.
# }
#
# \details{
#  Since the \pkg{affy} package is making use of special CDF environment
#  packages, this method will warn if the needed package is missing and
#  explain that \pkg{affy} will later try to download and install it
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
setMethodS3("extractAffyBatch", "AffymetrixCelSet", function(this, ..., verbose=FALSE) {
  require(affy) || throw("Package not loaded: affy");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # 
  cdf <- getCdf(this);
  cdfPkgName <- cleancdfname(getChipType(cdf));
  if (!suppressWarnings(require(cdfPkgName, character.only=TRUE))) {
    warning("CDF enviroment package '", cdfPkgName, "' not installed. The 'affy' package will later try to download from Bioconductor and install it.");
  }

  filenames <- getPathnames(this);
  verbose && enter(verbose, "Creating AffyBatch from ", length(filenames), " CEL files");
  verbose && cat(verbose, "Filenames: ", paste(filenames, collapse=", "));
  sampleNames <- getNames(this);
  verbose && cat(verbose, "Sample names: ", paste(sampleNames, collapse=", "));

  res <- ReadAffy(filenames=filenames, sampleNames=sampleNames, ..., verbose=as.logical(verbose));

  verbose && exit(verbose);
  res;
}) # extractAffyBatch()


############################################################################
# HISTORY:
# 2006-10-02
# o Created. A first small step toward an interface to Bioconductor.
############################################################################
