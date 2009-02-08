setMethodS3("fitRoc", "SegmentedCopyNumbers", function(this, states=NULL, recall=states[2], ...) {
  # Extract by state?
  if (!is.null(states)) {
    this <- extractSubsetByState(this, states=states);
  }

  # Extract data
  data <- as.data.frame(this);
  rm(this);

  # Recall?
  if (!is.null(recall)) {
    data$state <- as.integer(is.element(data$state, recall));
  }

  fitRoc(truth=data$state, data=data$cn, ...);
})


setMethodS3("findRocTpAtFp", "SegmentedCopyNumbers", function(this, states=NULL, recall=states[2], ...) {
  # Extract by state?
  if (!is.null(states)) {
    this <- extractSubsetByState(this, states=states);
  }

  # Extract data
  data <- as.data.frame(this);
  rm(this);

  # Recall?
  if (!is.null(recall)) {
    data$state <- as.integer(is.element(data$state, recall));
  }

  findRocTpAtFp(truth=data$state, data=data$cn, ...);
})



############################################################################
# HISTORY:
# 2009-02-07
# o Created.
############################################################################
