source("init.R")

verbose <- Arguments$getVerbose(TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("ds")) {
  path <- "chip_data3/Xba/";
  ds <- AffymetrixCelSet$fromFiles(path);
}

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Create a set of models to work with
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
if (!exists("models", mode="list")) {
  mergeStrands <- TRUE;
  models <- list(
    rma = RmaSnpPlm(ds, mergeStrands=mergeStrands),
    mbei = MbeiSnpPlm(ds, mergeStrands=mergeStrands),
    affine = AffineSnpPlm(ds, background=FALSE, mergeStrands=mergeStrands)
  )
}

units <- 55+1:300;

lapply(models, function(model) {
  print(model);
  t <- system.time({
    uu <- fit(model, units=NULL, force=FALSE, verbose=TRUE);
  })
  if (length(uu) > 0)
    print(t / length(uu));
})

log <- TRUE;
opar <- par(ask=(length(units) > 1));
for (unit in units) {
  for (kk in seq(along=models)) {
    model <- models[[kk]];
    ces <- getChipEffects(model);
    ceUnit <- ces[unit];
    y <- ceUnit[[1]];
    snpName <- names(ceUnit)[1];
    yA <- y[[1]]$theta;
    yB <- y[[2]]$theta;
    if (log) {
      yA <- log(yA, base=2);
      yB <- log(yB, base=2);
      lim <- c(6, 16);
    } else {
      lim <- c(0, 2^15);
    }
    if (kk == 1) {
      plot(NA, xlim=lim, ylim=lim, xlab="A", ylab="B", main=snpName);
      abline(a=0, b=1, lty=2);
    }
    points(yA, yB, col=kk, pch=19);
  }
}
par(opar);