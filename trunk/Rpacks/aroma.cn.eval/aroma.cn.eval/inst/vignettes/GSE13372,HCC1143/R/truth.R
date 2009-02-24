############################################################################
# Known CN aberrant regions
############################################################################
regions <- c(
  "GSM337641:Chr1@98-112",
  "GSM337641:Chr3@80-92",
  "GSM337641:Chr4@60.5-68.5",
  "GSM337641:Chr10@61-69",
  "GSM337641:Chr11@78.1-83",
  "GSM337641:Chr12@57-63"
);


truth <- function(x, chromosome, name, ...) {
  name <- gsub(",.*", "", name);
  res <- integer(length(x));
  state <- -1L;
  cps <- NULL;
  if (chromosome == 1) {
    cps <- c(103.8,107.7)*1e6;
    dx <- 200e3;
    state <- -1L;
  } else if (chromosome == 3) {
    cps <- c(85.3,91)*1e6;
    dx <- 100e3;
    state <- +1L;
  } else if (chromosome == 4) {
    cps <- c(63.4,65.8)*1e6;
    dx <- 50e3;
    state <- +1L;
  } else if (chromosome == 10) {
    cps <- c(60,65.3)*1e6;
    dx <- 100e3;
    state <- +1L;
  } else if (chromosome == 11) {
    cps <- c(78.1,80.2)*1e6;
    dx <- 100e3;
    state <- +1L;
  } else if (chromosome == 12) {
    cps <- c(56.0,59.9)*1e6;
    dx <- 100e3;
    state <- +1L;
  }

  if (length(cps) > 0) {
    res[cps[1] <= x & x < cps[2]] <- state;
    dx <- rep(dx, length.out=2);
    for (kk in seq(along=cps)) {
      res[cps[kk]-dx[kk] <= x & x < cps[kk]+dx[kk]] <- NA;
    }
  }
  res;
} # truth()



############################################################################
# HISTORY:
# 2009-02-23
# o Created.
############################################################################
