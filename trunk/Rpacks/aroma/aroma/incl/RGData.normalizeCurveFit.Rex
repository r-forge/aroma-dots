# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Main example code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# One array was scanned four times at four different PMT settings.
rg <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
setLayout(rg, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg, c("500V","600V","700V","800V"));
# Reorder the slides in scan order
keepSlides(rg, c("800V","500V","600V","700V"));
# Since scan "500V" seems to be an outlier, exclude it!
removeSlides(rg, slide="500V")

rg0 <- clone(rg);

# Calibrate using x <- (y-a)/b. A better alternative is to projecting the 
# fitted line onto (1,1,...,1), but in this example we want to show that
# the scan after translation and rescaling are very similar, but with noise.
fit <- calibrateMultiscan(rg, project=FALSE);

avg <- median
R <- apply(rg$R, MARGIN=1, FUN=avg, na.rm=TRUE)
G <- apply(rg$G, MARGIN=1, FUN=avg, na.rm=TRUE)
R <- R / 8; G <- G/8;
rg <- RGData(R=R,G=G,layout=getLayout(rg))
rm(R,G)

# Since we used project=TRUE all scans are identical now.
keepSlides(rg, slide=1)

Device$set(2, height="108%");
subplots(9, nrow=3);

# Plot such that there is a right angle between the log(R) and log(G) axes.
Alim <- c(-2,18)
Mlim <- c(-1,1)*abs(diff(Alim))

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main="multiscan calibrated");

# Curve-fit using splines is much faster than using loess and give
# basically the same estimate.
fit <- normalizeCurveFit(rg, method="spline")

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main="calibrated + normalized");
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main="calibrated + normalized");


rg <- clone(rg0)
keepSlides(rg, slide="600V")

ma <- as.MAData(rg);
pmt <- paste("PMT", getSlideNames(ma));
plot(ma, xlim=Alim, ylim=Mlim, main=pmt);

cat("Signals (", pmt, ") before normalization:\n", sep="");
print(summary(rg))
fit <- normalizeCurveFit(rg, method="spline")
cat("Signals (", pmt, ") after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)

ma <- as.MAData(rg);
plot(ma, xlim=Alim, ylim=Mlim, main=paste(pmt, "normalized"));
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main=paste(pmt, "normalized"));


rg <- clone(rg0)

ma <- as.MAData(rg);
plot(NA, xlim=Alim, ylim=Mlim, main="All scans");
for (kk in seq(ma))
  points(ma, slide=kk, col=kk+1, pch=".");

cat("Signals before normalization:\n", sep="");
print(summary(rg))
fit <- normalizeCurveFit(rg, method="spline")
cat("Signals after normalization:\n", sep="");
print(summary(rg))
cat("based on the estimates:\n");
str(fit)

# To give evidence that all signals now contains the same bias 
# and that the relative scale between all scans and all channels
# is one we rescale all signals by the same factor. The only affect
# this will have on the M vs A plot is that the data is shifted
# along the A dimension. Note that there is no rescaling the 
# log-ratios!
#scale <- 1/2^5;   # Shift such that A <- A - 5;
#for (ch in c("R", "G"))
#  rg[[ch]] <- scale*rg[[ch]];

ma <- as.MAData(rg);
plot(NA, xlim=Alim, ylim=Mlim, main="All scans normalized");
for (kk in seq(ma))
  points(ma, slide=kk, col=kk+1, pch=".");
drawRGGrid(maxSignal=Alim[2])

plotPairDensities(rg, ylim=c(0,0.4), main="All scans normalized");

