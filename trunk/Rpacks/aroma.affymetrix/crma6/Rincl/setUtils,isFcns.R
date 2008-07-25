# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# isNnn() functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
isSmoothedSet <- function(sets, ...) {
  (regexpr("h=", names(sets)) != -1);
}

isRawSet <- function(sets, ...) {
  (regexpr("h=", names(sets)) == -1);
}

isSubsetSet <- function(sets, subset=NULL, ...) {
  if (identical(subset, "all")) {
    pattern <- paste("subset=", sep="");
    (regexpr(pattern, names(sets)) == -1);
  } else {
    pattern <- paste("subset=", subset, sep="");
    (regexpr(pattern, names(sets)) != -1);
  }
}

isRawDataSet <- function(sets, ...) {
  (regexpr(",", names(sets)) == -1);
}
