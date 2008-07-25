setConstructorS3("SmartList", function(...) {
  this <- list(...);
  class(this) <- c("SmartList", class(this));
  this;
})

setMethodS3("[[<-", "SmartList", function(this, name, value) {
  oldValue <- this[[name]];
  if (identical(value, oldValue)) {
#    cat("Nothing changed!\n");
    invisible(this);
  } else {
    NextMethod("[[<-");
  }
})

setMethodS3("$<-", "SmartList", function(this, name, value) {
#  cat("$<-....\n");
  UseMethod("[[<-");
})

setMethodS3("[", "SmartList", function(this, ...) {
  class <- class(this);
  this <- NextMethod("[");
  class(this) <- class;
  this;
})


############################################################################
# HISTORY:
# 2008-07-24
# o Created in order to avoid a copy of a list if an assignment makes no
#   difference.
############################################################################
