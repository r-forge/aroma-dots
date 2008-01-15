library(aroma)
source("RGData.R")

if (!exists("rg")) {
  # Create a raw data object from the preexisting example data in
  # the sma package.
  SMA$loadData("mouse.data")
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  
  # Get the signal (here by default non-background corrected)
  ma <- getSignal(raw)
  
  # Transform (M,A) into (R,G)
  rg <- as.RGData(ma)
}

printtips <- getPrinttipGroups(layout);
spots <- getSpots(printtips);

slide <- 1;
R <- rg$R[,slide];
G <- rg$G[,slide];

incl <- unlist(spots[1:11]);
R <- R[incl]; G <- G[incl];

ok <- is.finite(R) & is.finite(G);
R <- R[ok];
G <- G[ok];


require("multidim") || throw("Package 'multidim' not found.");
require("MASS") || throw("Package 'MASS' not found.");

maxIterAboveBelow <- 20;
maxIter <- 200;
reps <- 0.01;

xy <- cbind(R,G);

lim <- c(0,1e4);
plot(xy, pch=".", xlim=lim, ylim=lim);

w1 <- rep(1, length=nrow(xy));
wCenter <- w1;
for (iter in 1:maxIter) {
  # center line
  fitCenter <- acp(xy, wt=wCenter, reduc=FALSE);
  iCenter <- fit$cmpr[,1];
  rCenter <- fit$cmpr[,2];
  A <- fitCenter$vectors;
  x0 <- fitCenter$moy;
  m <- x0[2]-A[2,1]/A[1,1]*x0[1];
  k <- A[2,1]/A[1,1];

  abline(a=m,b=k, col="red");

  # above line
  above <- (rCenter > 0);
# points(xyAbove, pch=".", col="red");
  xyAbove <- xy[above,];
  wAbove <- rep(1, length=nrow(xyAbove));
  for (iterAbove in 1:maxIterAboveBelow) {
    fitAbove <- rlm(xyAbove, wt=wAbove, reduc=FALSE);
    wAbove <- 1/abs(fitAbove$cmpr[,2]+reps);
    A <- fitAbove$vectors;
    x0 <- fitAbove$moy;
    mA <- x0[2]-A[2,1]/A[1,1]*x0[1];
    kA <- A[2,1]/A[1,1];
  }

  abline(a=mA,b=kA, col="red", lty=2);

  # below line
  xyBelow <- xy[!above,];
# points(xyBelow, pch=".", col="green");
  wBelow <- rep(1, length=nrow(xyBelow));
  for (iterAbove in 1:maxIterAboveBelow) {
    fitBelow <- acp(xyBelow, wt=wBelow, reduc=FALSE);
    wBelow <- 1/abs(fitBelow$cmpr[,2]+reps);
    A <- fitBelow$vectors;
    x0 <- fitBelow$moy;
    mB <- x0[2]-A[2,1]/A[1,1]*x0[1];
    kB <- A[2,1]/A[1,1];
  }

  abline(a=mB,b=kB, col="green", lty=2);

  # Scale function
  # Rescaled residuals
  # The weight for the center line should be
  wCenter <- 1/abs(rCenter+reps);

  cat(sprintf("%4d. m=%.3f, k=%.3f, ", iter, m, k));
  cat(sprintf("mA=%.3f, kA=%.3f, ", mA, kA));
  cat(sprintf("mB=%.3f, kB=%.3f\n", mB, kB));
}


abline(a=m,b=k, col="black", lty=2, lwd=2);
abline(a=mA,b=kA, col="black", lty=2, lwd=2);
abline(a=mB,b=kB, col="black", lty=2, lwd=2);




