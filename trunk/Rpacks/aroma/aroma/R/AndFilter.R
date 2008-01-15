#########################################################################/**
# @RdocClass AndFilter
#
# @title "The AndFilter class"
#
# \description{
#  @classhierarchy
#
#   An AndFilter is a ParallelFilter that passes through indices that are
#   outputs from \emph{all} input filters connected to this filter.
#   This corresponds to the logical operator \emph{AND}. The filter could
#   also have been called an \emph{intersection filter}.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{The input @see "Filter"s to be connected to.}
# }
#
# \section{Fields and Methods}{
#  @allmethods
# }
#
#
# @author
#
# \examples{\dontrun{See help(ParallelFilter) for an example.}}
#
# \seealso{
#   @see "OrFilter" and @see "ParallelFilter".
# }
#*/#########################################################################
setConstructorS3("AndFilter", function(...) {
  extend(ParallelFilter(...), "AndFilter"
  )
});




#########################################################################/**
# @RdocMethod getIndex
#
# @title "Gets indices accepted by this filter"
#
# @synopsis
#
# \description{
#  @get "title".
# }
#
# \value{
#   Returns a @vector of @integers.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# }
#*/######################################################################### 
setMethodS3("getIndex", "AndFilter", function(this) {
  filters <- this$filters;
  filter <- filters[[1]];
  idx <- getIndex(filter);
  maxIdx <- attr(idx, "max");
  for (filter in filters[-1]) {
    incl <- getIndex(filter);
    if (maxIdx < attr(incl, "max"))
      maxIdx <- attr(incl, "max");
    idx <- intersect(idx, incl);
  }
  attr(idx, "max") <- maxIdx;
  idx;
})



############################################################################
# HISTORY:
# 2003-05-06
# o Added Rdoc comments.
# o The AndFilter was mistakenly set to be abstract.
# 2002-02-26
# o Modified code to make use of setMethodS3's.
# o Cut'n'pasted from ParallelFilter.
############################################################################
