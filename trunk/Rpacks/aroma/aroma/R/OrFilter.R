#########################################################################/**
# @RdocClass OrFilter
#
# @title "The OrFilter class"
#
# \description{
#  @classhierarchy
#
#   An OrFilter is a ParallelFilter that passes through indices that are
#   outputs from \emph{any} input filters connected to this filter.
#   This corresponds to the logical operator \emph{OR}. The filter could
#   also have been called an \emph{union filter}.
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
#   @see "AndFilter" and @see "ParallelFilter".
# }
#*/#########################################################################
setConstructorS3("OrFilter", function(...) {
  extend(ParallelFilter(...), "OrFilter"
  )
})


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
setMethodS3("getIndex", "OrFilter", function(this) {
  filters <- this$filters;
  filter <- filters[[1]];
  idx <- getIndex(filter);
  maxIdx <- attr(idx, "max");
  for (filter in filters[-1]) {
    incl <- getIndex(filter);
    if (maxIdx < attr(incl, "max"))
      maxIdx <- attr(incl, "max");
    idx <- union(idx, incl);
  }
  attr(idx, "max") <- maxIdx;
  idx;
})



############################################################################
# HISTORY:
# 2003-06-23
# o BUG FIX: RdocClass was AndFilter.
# 2003-05-06
# o Added Rdoc comments.
# o The OrFilter was mistakenly set to be abstract.
# 2002-02-26
# o Modified code to make use of setMethodS3's.
# o Cut'n'pasted from ParallelFilter.
############################################################################
