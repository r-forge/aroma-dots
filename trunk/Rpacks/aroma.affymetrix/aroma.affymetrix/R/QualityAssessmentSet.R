###########################################################################/**
# @RdocClass QualityAssessmentSet
#
# @title "The QualityAssessmentSet class"
#
# \description{
#  @classhierarchy
#
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "QualityAssessmentFile":s.}
#   \item{tags}{A @character @vector of tags.}
#   \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{
#   Ken Simpson (ksimpson[at]wehi.edu.au).
# }
#*/###########################################################################
setConstructorS3("QualityAssessmentSet", function(files=NULL, tags="*", ...) {
  
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));

    # Update default tags
    tags[tags == "*"] <- "QA";
  }

  extend(AffymetrixCelSet(files=files, ...), "QualityAssessmentSet");
  
})



#setMethodS3("fromFiles", "QualityAssessmentSet", function(static, path, pattern="[.]calls$", ..., fileClass="GenotypeCallFile", verbose=FALSE) {
#  # Argument 'verbose':
#  verbose <- Arguments$getVerbose(verbose);
#  if (verbose) {
#    pushState(verbose);
#    on.exit(popState(verbose));
#  }#
#
#  # S3 method dispatch does not work for static methods.
#  ds <- fromFiles.AffymetrixFileSet(static, path=path, pattern=pattern, ..., fileClass=fileClass, verbose=less(verbose));#
#
#  # Use the same CDF object for all CEL files.
#  setCdf(ds, getCdf(ds));#
#
#  ds;
#}, static=TRUE)

