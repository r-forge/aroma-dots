source("init.R")

verbose <- Arguments$getVerbose(TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Specify the dataset to be used
ds <- AffymetrixCelSet$fromFiles("chip_data/Hind/")
print(ds)
cdf <- getCdf(ds)

cel <- ds[[1]]

# Plot stdvs as a function of intensities for a subset of the cells
subset <- seq(from=1, to=nbrOfCells(cdf), length=1e4)
X <- getFields(cel, indices=subset, fields=c("intensities", "stdvs"))
X <- as.matrix(X)
xlim <- c(0,65535) / 16
ylim <- xlim
plot(X, xlim=xlim, ylim=ylim, xlab="y", ylab=expression(sigma), main=getName(cel))
usr <- par("usr")
abline(a=0, b=1, col="blue", lty=2)
points(0,0, pch="+", cex=2, col="blue")
dim <- paste(getDimension(cdf), collapse="x");
label <- sprintf("Chip type: %s", getChipType(cdf))
text(x=usr[2], y=usr[4], labels=label, adj=c(1,-0.5), cex=0.8, xpd=TRUE)

# Fit a line robustly
ok <- (X[,"intensities"] <= 2/3*(2^16-1))
X <- X[ok,]
fit <- fitIWPCA(X, constraint=0)  # From aroma.light
print(fit)
dx <- diff(par("usr")[1:2])
x <- fit$a[1] + c(0,1) * fit$b[1] * dx
y <- fit$a[2] + c(0,1) * fit$b[2] * dx
lines(x, y, col="red", lwd=2)
points(fit$a[1], fit$a[2], cex=2, lwd=2,  col="red")
