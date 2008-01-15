# Load some example data
SMA$loadData("mouse.data")
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
raw <- RawData(mouse.data, layout=layout)

# Get the signal
ma <- getSignal(raw, bgSubtract=TRUE)

subplots(4)

plotMvsM(ma, slides=c(1,2), xlim=c(-1,1)*3, pch=".")
points(ma, slides=c(3,4), col="red", pch=".")
points(ma, slides=c(5,6), col="blue", pch=".")
abline(a=0, b=1, lty=2)
title(main="Three different (M,M) pairs")


normalizeWithinSlide(ma, method="l")
normalizeAcrossSlides(ma)

plotMvsM(ma, slides=c(1,2), xlim=c(-1,1)*3, pch=".")
points(ma, slides=c(3,4), col="red", pch=".")
points(ma, slides=c(5,6), col="blue", pch=".")
abline(a=0, b=1, lty=2)
title(main="Same after lowess normalization")
