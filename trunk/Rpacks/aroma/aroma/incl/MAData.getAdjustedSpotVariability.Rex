SMA$loadData("mouse.data")
# Keep only slides of treatment 2
mouse.data <- lapply(mouse.data, FUN=function(x) x[,4:6])
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
raw <- RawData(mouse.data, layout=layout)

# Create four sets of slides where slide 2 and 4 are shifted R=G=a
a <- 2^11
ma <- list()
ma[[1]] <- getSignal(raw, bgSubtract=TRUE)
ma[[2]] <- clone(ma[[1]])
shift(ma[[2]], R=a, G=a)
ma[[3]] <- clone(ma[[1]])
normalizeWithinSlide(ma[[3]], method="p")
ma[[4]] <- clone(ma[[3]])
shift(ma[[4]], R=a, G=a)

# Calculates the (adjusted) spot variability 
d <- list()
for (k in 1:length(ma))
  d[[k]] <- getAdjustedSpotVariability(ma[[k]])

# Four plots
subplots(4)
# Plot 1 and 2: Non-normalized and print-tip normalized slides
# where each exists in a shifted and a non-shifted version
Alim <- Mlim <- NA
for (k in 1:length(ma)) {
  Alim <- range(c(Alim, ma[[k]]$A), na.rm=TRUE)
  Mlim <- range(c(Mlim, ma[[k]]$M), na.rm=TRUE)
}
plot(ma[[1]], xlim=Alim, ylim=Mlim, col=1)
points(ma[[2]], col=2)
plot(ma[[3]], xlim=Alim, ylim=Mlim, col=3)
points(ma[[4]], col=4)

# Plot 3: Densities of the (non-adjusted) spot variabilities
ymax <- xmax <- NA
ds <- list()
for (k in 1:length(ma)) {
  ds[[k]] <- density(na.omit(d[[k]]$d))
  ymax <- max(ymax, ds[[k]]$y, na.rm=TRUE)
}
for (k in 1:length(ma))
  xmax <- max(xmax, ds[[k]]$x[ds[[k]]$y > 0.01*ymax], na.rm=TRUE)
xlim <- c(0,xmax)
ylim <- c(0,ymax)
plot(NA, xlim=xlim, ylim=ylim, xlab="variability", ylab="density", 
					    main="Spot variabilities")
for (k in 1:length(ma))
  lines(ds[[k]], col=k)

# Plot 4: Densities of the *adjusted* spot variabilities
ymax <- xmax <- NA
ds <- list()
for (k in 1:length(ma)) {
  ds[[k]] <- density(na.omit(d[[k]]$dw))
  ymax <- max(ymax, ds[[k]]$y, na.rm=TRUE)
}
for (k in 1:length(ma))
  xmax <- max(xmax, ds[[k]]$x[ds[[k]]$y > 0.01*ymax], na.rm=TRUE)
xlim <- c(0,xmax)
ylim <- c(0,ymax)
plot(NA, xlim=xlim, ylim=ylim, xlab="variability", ylab="density", 
				   main="Adjusted spot variabilities")
for (k in 1:length(ma))
  lines(ds[[k]], col=k)                                                 
