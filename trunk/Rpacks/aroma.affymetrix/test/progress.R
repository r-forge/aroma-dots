library(aroma.affymetrix)
source("init.R")

pngDev <- System$findGraphicsDevice();

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
writeDoneBar <- function(done, filename="done.png", path="images", col=c(0,palette()[1]), width=1600) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path);

  if (writeImage) {
    nDone <- sum(done);
    progress <- nDone / length(done);
    if (0 < progress & progress < 1) {
      z <- matrix(as.integer(done), ncol=width);
      z <- colMeans(z);
      # Threshold such to achieve 'progress' saturation.
      qz <- quantile(z, probs=1-progress);
      z <- matrix(round((z >= qz)), nrow=1);
    } else {
      z <- matrix(progress, nrow=1, ncol=length(done));
    }

    # Write "done" bar
    pngDev(pathname, width=width, height=64, bg="white");
    on.exit(dev.off());
    par(mar=c(0,0,0,0));
    image270(z=z, zlim=c(0,1), col=col, axes=FALSE);
  }

  invisible(pathname);
} # writeDoneBar()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data-set specific data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!exists("ce")) {
  pathname <- file.choose();
  ce <- ChipEffectFile$fromFile(pathname);
  rm(todo);
}
print(ce);
cdf <- getCdf(ce);
print(cdf);
gi <- getGenomeInformation(cdf);
print(gi);
if (!exists("todo")) {
  todo <- findUnitsTodo(ce, verbose=TRUE);
  writeImage <- TRUE;
} else {
  writeImage <- FALSE;
}
str(todo);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Identify progress for whole genome and per chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
unitsDone <- list();

# Whole genome...
done <- rep(TRUE, nbrOfUnits(cdf));
done[todo] <- FALSE;
n <- length(done);
count <- sum(done);
unitsDone[["All"]] <- list(label="All", n=n, done=done, count=count, progress=count/n);

# Chromosome by chromosome...
allChromosomes <- c(1:22, "X");
for (chr in allChromosomes) { 
  label <- sprintf("%02d", match(chr, allChromosomes));
  label <- chr;
  print(label);
  units <- getUnitIndices(gi, chromosome=chr);
  n <- length(units);
  count <- sum(done[units]);
  unitsDone[[chr]] <- list(label=label, n=n, done=done[units], count=count, progress=count/n);
}

# Create done bars
for (kk in seq(along=unitsDone)) {
  set <- unitsDone[[kk]];
  print(set$label);
  pathname <- writeDoneBar(set$done, filename=sprintf("done%s.png", set$label));
  set$image <- pathname;
  unitsDone[[kk]] <- set;
}
