# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Number of loci
J <- 1000;

y <- rnorm(J, sd=1/2);
y[200:300] <- y[200:300] + 1;
y[650:800] <- y[650:800] - 1;
x <- sort(runif(length(y), max=length(y)));


cn <- RawCopyNumbers(y, x);
cn$fullname <- "Sample #1";
print(cn);

plot(cn, col="#aaaaaa");

cnS <- binnedSmoothing(cn, by=10);
print(cnS);
lines(cnS, col="black", lwd=3);

fit <- segmentByCBS(cn);
drawCnRegions(fit, col="red", lwd=3);

fit <- segmentByGLAD(cn);
drawCnRegions(fit, col="blue", lwd=3);

fit <- segmentByHaarSeg(cn);
drawCnRegions(fit, col="orange", lwd=3);
