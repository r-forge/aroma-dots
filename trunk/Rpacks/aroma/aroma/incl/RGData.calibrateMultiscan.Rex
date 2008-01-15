# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Main example code
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# One array was scanned four times at four different PMT settings.
rg0 <- RGData$read("PMT-RGData.dat", path=system.file("data-ex", package="aroma.light"));
setLayout(rg0, Layout(4,4,17,17));       # Not really necessary!
setSlideName(rg0, c("500V","600V","700V","800V"));
# Reorder the slides (to make the M vs A plot more clear!)
keepSlides(rg0, c("800V","500V","600V","700V"));

for (ch in c("R", "G")) {
  signalWeights <- SignalWeights$fromSaturatedSignals(rg0[[ch]]);
  setSignalWeights(rg0, channel=ch, weights=signalWeights);
  rm(signalWeights);
}

# Calibrate using x <- (y-a)/b. For illustration purposed, we want to 
# keep all scans so we do not calculate the average scan here. 
rgC <- clone(rg0);
fit <- calibrateMultiscan(rgC, average=NULL, project=FALSE);
str(fit);

Device$set(2, height="108%");
subplots(9);

Alim <- c(0,20);
Mlim <- c(-1,1)*abs(diff(range(Alim)));

# Same scale on calibrated and non-calibrated data
xlim <- c(4,18); ylim <- c(-1,1)*3.5;

for (ch in c("R", "G")) {
  plotPairwiseScans(rg0, channel=ch, xlim=Alim, ylim=Mlim, main="non-calibrated");
  drawRGGrid(afit=fit[[ch]]);
}

plotPairDensities(rg0, main="non-calibrated");

# It is clear from the calibrated plots that the 500V scan is special. 
# By visual inspection of the density plots the following manual 
# correction makes their densities "nicer". Uncomment to test.
# rg$R[,1] <- rg$R[,1] - 94; rg$G[,1] <- rg$G[,1] - 54;

for (ch in c("R", "G")) {
  plotPairwiseScans(rgC, channel=ch, xlim=Alim, ylim=Mlim, main="calibrated");
  drawRGGrid(maxSignal=Alim[2]);
}

plotPairDensities(rgC, main="calibrated");


# Now when we have shown the effect of multiscan calibration, we
# calibrate and calculate the median scan as done in practice.
rgCAvg <- clone(rg0);
calibrateMultiscan(rgCAvg);
setSlideNames(rgCAvg, "calib");


# Put both raw and average calibrated data into the same MAData object
rgN <- clone(rgCAvg);
fit <- normalizeAffine(rgN);

# Try also this to see the effect of lowess normalization
maN <- as.MAData(rgN); normalizeWithinSlide(maN, method="l"); rgN <- as.RGData(maN); 
rg <- rgN;
append(rg, rgCAvg);
setSlideNames(rg, c("normalized", "raw"));
rm(rgN,rgCAvg);

# Plot M vs A for raw data. From the overlaid grid (according to the 
# estimated affine parameters), one can see that the normalization
# method indeed allows some signals to fall out of the positive
# quadrant, that is, the data points outside the grid will become 
# negative values (and therefore not visible in the next (A,M)-plot).
ma <- as.MAData(rg);
Alim <- c(0,ceiling(max(ma$A,na.rm=TRUE))); 
Mlim <- c(-1,1)*abs(diff(range(Alim)));
plot(ma, slide="raw", xlim=Alim, ylim=Mlim);
drawRGGrid(afit=fit[[1]], drawCurve=TRUE);
plot(ma, slide="normalized", xlim=Alim, ylim=Mlim);
drawRGGrid(maxSignal=Alim[2]);

# Plot M vs M for normalized and raw data. Stratify on observations
# close to A = [0.0,0.5,1.0,...,16.0] to better see the intenstity 
# dependency of the log-ratio bias.
ma <- as.MAData(rg);

isCloseTo <- function(x, value=0, dx=0.05) (x > value-dx & x < value+dx);
incl <- rep(NA, nbrOfSpots(rg)); 
for (aa in seq(0,16, by=0.5))
  incl <- incl | isCloseTo(ma$A[,2], value=aa, dx=0.05);

# Shift the plot to maximize resolution, but keep aspect ratios.
Mlim <- c(-7,4);
plotMvsM(ma, slides=1:2, pch=176, cex=0.5, include=incl, xlim=Mlim); 
abline(a=0,b=1, col="lightblue");

