###########################################################################/**
# @RdocClass AptBrlmm
#
# @title "The AptBrlmm class"
#
# \description{
#  @classhierarchy
#
#  This class represents an interface to the BRLMM algorithm implemented in
#  the Affymetrix Power Tools (APT) software [1].
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \details{
#  Note that the methods of this class are calling the APT [1] binaries,
#  which must be installed on the system.
# }
#
# @author
#
# \seealso{
# }
#
# \references{
#  [1] Affymetrix, Affymetrix Power Tools (APT) software, Dec 2006.
#      \url{http://www.affymetrix.com/support/developer/powertools/index.affx}
#      \cr
# }
#*/###########################################################################
setConstructorS3("AptBrlmm", function(...) {
  extend(Object(), "AptBrlmm",
    ...
  )
})


setMethodS3("as.character", "AptBrlmm", function(this, ...) {
  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("RAM: %.2fMb", objectSize(this)/1024^2));
  class(s) <- "GenericSummary";
  s;
})


###########################################################################/**
# @RdocMethod getRootPath
#
# @title "Gets the root path of this model"
#
# \description{
#  @get "title", which is the class name of the model preceeded by the
#  string \code{"model"}.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getRootPath", "AptBrlmm", function(this, ...) {
  # Default root path
  paste("model", class(this)[1], sep="");
})



###########################################################################/**
# @RdocMethod getName
#
# @title "Gets the name of the output summarized data set"
#
# \description{
#  @get "title", which is the same as the input data set.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getName", "AptBrlmm", function(this, ...) {
  stop("TO DO");
})


###########################################################################/**
# @RdocMethod getTags
#
# @title "Gets the tags of the output data set"
#
# \description{
#  @get "title", which equals the tags of the input data set plus the tags
#  of this model.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character @vector.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getTags", "AptBrlmm", function(this, ...) {
  stop("TO DO");
})



###########################################################################/**
# @RdocMethod getFullName
#
# @title "Gets the full name of the output data set"
#
# \description{
#  @get "title", which is the name with comma separated tags.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getFullName", "AptBrlmm", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})



###########################################################################/**
# @RdocMethod getPath
#
# @title "Gets the path of this model"
#
# \description{
#  @get "title" where the parameter files are stored.
#  If non-existing, then the directory is created.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#  Returns a @character string.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getPath", "AptBrlmm", function(this, ...) {
  # Create the (sub-)directory tree for the dataset

  # Root path
  rootPath <- getRootPath(this);
  mkdirs(rootPath);

  # Full name
  fullname <- getFullName(this);

  # Chip type    
  ces <- getChipEffects(this);
  cdf <- getCdf(ces);
  chipType <- getChipType(cdf);
  chipType <- gsub("-monocell$", "", chipType);

  # The full path
  path <- filePath(rootPath, fullname, chipType, expandLinks="any");
  if (!isDirectory(path)) {
    mkdirs(path);
    if (!isDirectory(path))
      throw("Failed to create output directory: ", path);
  }

  path;
})


###########################################################################/**
# @RdocMethod getChipEffects
#
# @title "Gets the chip effects for this model"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an @see "SnpChipEffectSet".
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("getChipEffects", "AptBrlmm", function(this, ...) {
  stop("TO DO");
})


###########################################################################/**
# @RdocMethod getCdf
#
# @title "Gets the CDF structure for this model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
# }
#
# \value{
#  Returns an @see "AffymetrixCdfFile" object.
# }
#
# @author
#
# \seealso{
#   @seemethod "setCdf".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("getCdf", "AptBrlmm", function(this, ...) {
  getCdf(this$ces);
})


###########################################################################/**
# @RdocMethod fit
#
# @title "Fits the model"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{recalibrate}{If @TRUE, a second round of the fitting is done.}
#   \item{...}{Additional arguments.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns a @see "oligo::SnpCallSetPlus-class" object.
# }
#
# \details{
#   Simple benchmarking on shadowfax: 
#   For 90 CEPH Xba chips it takes ~90 minutes, i.e. 60 s/array.
# }
#
# @author
#
# \seealso{
#   @see "oligo::crlmm".
#   @seeclass
# }
#
# @keyword IO
#*/###########################################################################
setMethodS3("fit", "AptBrlmm", function(this, recalibrate=TRUE, ..., verbose=FALSE) {
  stop("TO DO");
})


############################################################################
# HISTORY:
# 2006-12-01
# o Created a skeleton.
############################################################################
