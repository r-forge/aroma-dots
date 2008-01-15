SMA$loadData("mouse.data")

layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
replicates <- GeneGroups(layout)
raw <- RawData(mouse.data, layout=layout)

# Get the signal (here by default non-background corrected)
ma <- getSignal(raw)
normalizeWithinSlide(ma, "m")
normalizeAcrossSlides(ma)

print(replicates)

cat("Calculating the standardized residuals...")
idx <- getSpots(replicates)
M <- ma$M     # To improve speed!
residuals <- lapply(idx, FUN=function(i) {
  x <- M[i,]; ok <- !is.na(x)
  c <- median(x[ok]); s <- mad(x[ok], center=c)
  x <- (x-c)/s
  attr(x, "center") <- c; attr(x, "scale") <- s
  x
})
for (k in seq(idx)) M[idx[[k]],] <- residuals[[k]]
cat("done\n")

ma.res <- clone(ma)
ma.res$M <- M
subplots(ma)
for (s in seq(ma))
  plotSpatial(ma.res, slide=s)


