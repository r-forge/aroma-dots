if (interactive()) savehistory();
library("aroma.core");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

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
title(main="Non-weighted segmentation");

cols <- c(CBS="red", GLAD="blue", HaarSeg="orange");

fit <- segmentByCBS(cn);
drawCnRegions(fit, col=cols["CBS"], lwd=3);

fit <- segmentByGLAD(cn);
drawCnRegions(fit, col=cols["GLAD"], lwd=3);

fit <- segmentByHaarSeg(cn);
drawCnRegions(fit, col=cols["HaarSeg"], lwd=3);

legend("topleft", horiz=TRUE, bty="n", lwd=3, col=cols, names(cols));
