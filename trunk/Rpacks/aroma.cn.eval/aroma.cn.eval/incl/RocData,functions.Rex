# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulate copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 100

# Number of channels
I <- 4

# Copy-neutral data
Y <- matrix(rnorm(I*J, sd=1/2), nrow=J, ncol=I)
T <- function(x, channel, ...) {
  t <- integer(length(x));
  if (is.element(channel, c(2,3))) {
    t[30 < x & x <= 50] <- +1;
    t[(25 < x & x <= 30) | (50 < x & x <= 55)] <- NA;
  }
  if (is.element(channel, c(1,3))) {
    t[65 < x & x <= 90] <- -1;
  }
  t;
} # T()

rd <- CopyNumberRocData(data=Y, truth=T)
Y <- getData(rd);
isT <- hasState(rd, +1);
Y[isT] <- Y[isT] + 3
isT <- hasState(rd, -1);
Y[isT] <- Y[isT] - 3
setData(rd, Y);

devSet(2);
plotTracks(rd);


devSet(3);
subset <- seq(from=1, to=nbrOfLoci(rd), by=3);
rd2 <- extractSubset(rd, rows=subset);
plotTracks(rd2);

