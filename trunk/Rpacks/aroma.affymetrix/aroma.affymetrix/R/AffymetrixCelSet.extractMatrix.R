###########################################################################/**
# @set "class=AffymetrixCelSet"
# @RdocMethod extractMatrix
#
# @title "Extract data as a matrix for a set of arrays"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cells}{(The subset of cells to be matched.
#     If @NULL, all cells are considered.}
#   \item{...}{Not used.}
#   \item{field}{The field to be extracted.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#  Returns an JxK @double @matrix where J is the number of units, 
#  and K is the number of arrays.
#  The names of the columns are the names of the arrays.
#  No names are set for the rows.
#  The rows are ordered according to \code{cells} argument.
# }
#
# @author
#
# \seealso{
#   @seeclass
# }
#*/###########################################################################
setMethodS3("extractMatrix", "AffymetrixCelSet", function(this, cells=NULL, ..., field=c("intensities", "stdvs", "pixels"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  cdf <- getCdf(this);
  if (is.null(cells)) {
    ncells <- nbrOfCells(cdf);
  } else {
    cells <- Arguments$getIndices(cells, range=c(1,nbrOfCells(cdf)));
    ncells <- length(cells);
  }

  # Argument 'field':
  field <- match.arg(field);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Getting data for the array set");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Allocate return array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Allocating matrix");
  arrayNames <- getNames(this);
  nbrOfArrays <- length(arrayNames);
  df <- matrix(NA, nrow=ncells, ncol=nbrOfArrays);
  colnames(df) <- arrayNames;
  verbose && str(verbose, df);
  verbose && printf(verbose, "RAM: %.2fMB\n", object.size(df)/1024^2);
  verbose && exit(verbose);

  verbose && enter(verbose, "Optimize reading order");
  o <- order(cells);
  cells <- cells[o];
  o <- order(o);
  verbose && exit(verbose);
  
  # Garbage collect
  gc <- gc();
  verbose && print(verbose, gc);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get thetas from the samples
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving data");
  for (aa in seq_len(nbrOfArrays)) {
    verbose && printf(verbose, "Array %d,\n", aa);
    cf <- getFile(this, aa);
    df[o,aa] <- getData(cf, indices=cells, fields=field, 
                                           verbose=less(verbose))[[field]];
    if (aa %% 10 == 0) {
      # Garbage collect
      gc <- gc();
      verbose && print(verbose, gc);
    }
  } # for (aa in ...)
  verbose && exit(verbose);

  verbose && exit(verbose);

  df;
}) # extractMatrix()


############################################################################
# HISTORY:
# 2007-03-29
# o Created from ChipEffectSet.extractMatrix.R.
############################################################################
