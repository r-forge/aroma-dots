# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# True CN states
stateFcn <- function(x, ...) {
  states <- integer(length(x))
  states[200 <=x & x <= 500] <- -1L
  states[600 <=x & x <= 900] <- +1L
  states
}

# Number of loci
J <- 1000

y <- rnorm(J, sd=1/2)
x <- 1:length(y)
for (state in c(-1,+1)) {
  idxs <- (stateFcn(x) == state)
  y[idxs] <- y[idxs] + 0.5*state
}


cnR <- SegmentedCopyNumbers(y, x, states=stateFcn)
print(cnR)

states <- c(0,1)
cn <- extractSubsetByState(cnR, states=states)
print(cn)
plot(cn, ylim=c(-3,3))


fit <- fitRoc(cn, recall=states[1])
plot(fit$roc, type="l", lwd=4, col=1)
abline(a=0, b=1, lty=3)

bys <- seq(from=1, to=10, by=3)
print(bys)
for (kk in seq(along=bys)) {
  by <- bys[kk]
  cnS <- binnedSmoothingByState(cn, by=by)
  cnS <- extractSubsetByState(cnS, states=states)
  print(table(getStates(cnS)))
  fit <- fitRoc(cnS, recall=states[1])
  lines(fit$roc, lwd=4, col=kk)
}
