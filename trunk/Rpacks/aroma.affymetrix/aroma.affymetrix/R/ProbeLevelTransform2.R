###########################################################################/**
# @RdocClass ProbeLevelTransform2
#
# @title "The ProbeLevelTransform2 class"
#
# \description{
#  @classhierarchy
#
#  This abstract class is specialized from @see "ProbeLevelTransform" and
#  provides methods to identify subsets and types of probes that are used
#  for fitting and/or updating the signals.
# }
# 
# @synopsis 
#
# \arguments{
#   \item{dataSet}{A @see "AffymetrixCelSet".}
#   \item{...}{Arguments passed to the constructor of 
#     @see "ProbeLevelTransform".}
#   \item{subsetToUpdate}{The probes to be updated.
#     If @NULL, all probes are considered.}
#   \item{typesToUpdate}{Types of probes to be updated.}
#   \item{typesToFit}{Types of probes to be used when fitting the model.}
#   \item{subsetToFit}{The units from which the normalization curve should
#     be estimated.  If @NULL, all are considered.}
#   \item{shift}{An optional amount to shift data before fitting and updating.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"  
# }
# 
# @author
#*/###########################################################################
setConstructorS3("ProbeLevelTransform2", function(dataSet=NULL, ..., typesToUpdate="pm", subsetToUpdate=NULL, typesToFit=typesToUpdate, subsetToFit="-XY", shift=0) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  extraTags <- NULL;

  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    if (!inherits(dataSet, "AffymetrixCelSet")) {
      throw("Argument 'dataSet' is not an AffymetrixCelSet object: ", 
                                                          class(dataSet)[1]);
    }

    # Argument 'typesToUpdate':
    if (!is.null(typesToUpdate)) {
      typesToUpdate <- match.arg(typesToUpdate, choices=c("pm", "mm", "pmmm"));
    }

    # Argument 'subsetToUpdate':
    if (is.null(subsetToUpdate)) {
    } else if (is.character(subsetToUpdate)) {
      if (subsetToUpdate %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToUpdate': ", subsetToUpdate);
      }
      extraTags <- c(extraTags, subsetToUpdate=subsetToUpdate);
    } else {
      df <- getFile(dataSet, 1);
      nbrOfCells <- nbrOfCells(df);
      subsetToUpdate <- Arguments$getIndices(subsetToUpdate, range=c(1, nbrOfCells));
      subsetToUpdate <- unique(subsetToUpdate);
      subsetToUpdate <- sort(subsetToUpdate);
    }

    # Argument 'typesToFit':
    if (!is.null(typesToFit)) {
      typesToFit <- match.arg(typesToFit, choices=c("pm", "mm", "pmmm"));
    }

    # Argument 'subsetToFit':
    if (is.null(subsetToFit)) {
    } else if (is.character(subsetToFit)) {
      if (subsetToFit %in% c("-X", "-Y", "-XY")) {
      } else {
        throw("Unknown value of argument 'subsetToFit': ", subsetToFit);
      }
      extraTags <- c(extraTags, subsetToFit=subsetToFit);
    } else {
      df <- getFile(dataSet, 1);
      nbrOfCells <- nbrOfCells(df);
      subsetToFit <- Arguments$getIndices(subsetToFit, range=c(1, nbrOfCells));
      subsetToFit <- unique(subsetToFit);
      subsetToFit <- sort(subsetToFit);
    }
  }

  # Argument 'shift':
  shift <- Arguments$getDouble(shift, disallow=c("NA", "NaN", "Inf")); 


  extend(ProbeLevelTransform(dataSet=dataSet, ...), "ProbeLevelTransform2",
    shift = shift,
    .typesToUpdate = typesToUpdate,
    .subsetToUpdate = subsetToUpdate,
    .typesToFit = typesToFit,
    .subsetToFit = subsetToFit,
    .extraTags = extraTags
  )
})


setMethodS3("clearCache", "ProbeLevelTransform2", function(this, ...) {
  # Clear all cached values.
  for (ff in c(".subsetToUpdateExpanded", ".subsetToFitExpanded")) {
    this[[ff]] <- NULL;
  }

  # Then for this object 
  NextMethod("clearCache", object=this, ...);
})


setMethodS3("getAsteriskTags", "ProbeLevelTransform2", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags", this, collapse=collapse, ...);

  # Extra tags?
  tags <- c(tags, this$.extraTags);

  # Add class-specific tags
  shift <- as.integer(round(this$shift));
  if (shift != 0) {
    tags <- c(tags, sprintf("%+d", shift));
  } 

  # Collapse?
  tags <- paste(tags, collapse=collapse);

  tags;
}, private=TRUE)




setMethodS3("getSubsetTo", "ProbeLevelTransform2", function(this, what=c("fit", "update"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what);

  field <- sprintf(".subsetTo%s", capitalize(what));
  fieldExpanded <- paste(field, "Expanded", sep="");
  typesField <- sprintf(".typesTo%s", capitalize(what));

  subset <- this[[field]];
  stratifyBy <- this[[typesField]];

  # Expand?
  if (is.character(subset)) {
    cells <- this[[fieldExpanded]];
    if (is.null(cells)) {
      dataSet <- getInputDataSet(this);
      cdf <- getCdf(dataSet);
      units <- subset;
      cells <- getSubsetOfCellIndices(cdf, units=units, stratifyBy=stratifyBy, ...);
      this[[fieldExpanded]] <- cells;
    }
  } else if (is.numeric(subset)) {
    cells <- subset;
  } else if (is.null(subset)) {
    dataSet <- getInputDataSet(this);
    cdf <- getCdf(dataSet);
    if (is.null(stratifyBy)) {
      cells <- seq(length=nbrOfCells(cdf));
    } else {
      cells <- getSubsetOfCellIndices(cdf, stratifyBy=stratifyBy, ...);
    }
    this[[fieldExpanded]] <- cells;
  } else {
    throw("Internal error. This statment should never be reached: ", mode(subset));
  }

  cells;
}, private=TRUE);


setMethodS3("getSubsetToUpdate", "ProbeLevelTransform2", function(this, ...) {
  getSubsetTo(this, what="update", ...);
}, protected=TRUE);


setMethodS3("getSubsetToFit", "ProbeLevelTransform2", function(this, ...) {
  getSubsetTo(this, what="fit", ...);
}, protected=TRUE);



setMethodS3("getParameters", "ProbeLevelTransform2", function(this, expand=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Get parameters from super class
  params <- NextMethod(generic="getParameters", object=this, expand=expand, ...);

  params <- c(params, list(
    typesToUpdate = this$.typesToUpdate,
    subsetToUpdate = this$.subsetToUpdate,
    typesToFit = this$.typesToFit,
    subsetToFit = this$.subsetToFit,
    shift = this$shift
  ));

  # Expand?
  if (expand) {
    verbose && enter(verbose, "Expanding subsets to fit and update");
    params$subsetToUpdate <- getSubsetToUpdate(this, verbose=less(verbose, 1));
    params$subsetToFit <- getSubsetToFit(this, verbose=less(verbose, 1));
    verbose && exit(verbose);
  }

  params;
}, private=TRUE)



############################################################################
# HISTORY:
# 2008-07-20
# o Extracted from BaseCountNormalization.R.
############################################################################
