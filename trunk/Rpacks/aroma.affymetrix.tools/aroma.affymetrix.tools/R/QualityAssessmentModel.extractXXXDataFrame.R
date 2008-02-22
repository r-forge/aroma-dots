setMethodS3("extractRleDataFrame","QualityAssessmentModel",
function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose && enter(verbose, "Getting data for the array set");
  data <- extractRle(this, ..., returnUgcMap=TRUE,verbose=less(verbose, 1));
  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;
  ugcMap <- as.data.frame(ugcMap);
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && enter(verbose, "Merging UGC map and extracted data");
  data <- cbind(ugcMap, data);
  verbose && exit(verbose);
  verbose && exit(verbose);
  data;
}
setMethodS3("extractNuseDataFrame","QualityAssessmentModel",
function(this, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose && enter(verbose, "Getting data for the array set");
  data <- extractNuse(this, ..., returnUgcMap=TRUE,verbose=less(verbose, 1));
  ugcMap <- attr(data, "unitGroupCellMap");
  attr(data, "unitGroupCellMap") <- NULL;
  ugcMap <- as.data.frame(ugcMap);
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);
  verbose && enter(verbose, "Merging UGC map and extracted data");
  data <- cbind(ugcMap, data);
  verbose && exit(verbose);
  verbose && exit(verbose);
  data;
}
