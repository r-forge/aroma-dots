log2neg <- function(x) { 
  x <- -x;
  x[x <= 0] <- NA; 
  log2(x); 
}

log2pos <- function(x) { 
  x[x <= 0] <- NA; 
  log2(x); 
}

log2abs <- function(x) { 
  x <- abs(x);
  x[x <= 0] <- NA; 
  log2(x); 
}

sqrtneg <- function(x) {
  x <- -x;
  x[x <= 0] <- NA; 
  sqrt(x); 
}

sqrtpos <- function(x) { 
  x[x <= 0] <- NA; 
  sqrt(x); 
}

sqrtabs <- function(x) { 
  x <- abs(x);
  x[x <= 0] <- NA; 
  sqrt(x); 
}


############################################################################
# HISTORY:
# 2007-02-14
# o Created.
############################################################################
