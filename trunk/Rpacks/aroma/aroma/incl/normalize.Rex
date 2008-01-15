library(R.io)   # File, FileOutputStream, HTMLPrintStream, ...
library(R.lang) # System

########################################################################
# 0. Crete a HTML report file
########################################################################
path <- File("C:/users/braju.com.R/public_html/com.braju.sma/work-spatial/CALLOW01")
mkdirs(path)
file <- File(path, "index.html")
fout <- FileOutputStream(file)
out <- HTMLPrintStream(fout)

writeDocType(out)
printTimestamp(out)
pushTag(out, "html")
pushTag(out, "head")
title <- "Normalization of the ApoA1 Knock-out data set";
printTag(out, "title", title)
popTag(out)
printTag(out, "meta", content="text/html; charset=iso-8859-1", "http-equiv"="Content-Type")
printTag(out, "meta", name="Author", content="[R] and R.classes, http://www.r-project.org")
pushTag(out, "body")
printComment(out, "B O D Y")

pushTag(out, "center")
printTag(out, "h1", title)
printTag(out, "hr")
popTag(out)

#source("MicroarrayData.R"); createClass("MicroarrayData")
#source("MAData.R"); createClass("MAData")

# Load the data
if (!exists("ma")) {
  data(MouseArray)
  
  layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
  raw <- RawData(mouse.data, layout=layout)
  
  # Get the signal
  ma <- getSignal(raw)
}

fig <- 1

########################################################################
# 1. Print order artifacts
########################################################################
printTag(out, "h1", "Print order artifacts")

imgfile <- File(path, "PrintorderArtifacts.png")
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)
  subplots(2, nrow=2)

  # If we order the data points in the same order as the spots where
  # printed on the cDNA microarray slide we see that both M and A
  # depend strongly on when they where printed.
  plotPrintorder(ma, slide=1)
  
  # Looking carefully at the data and using the knowledge that at each
  # turn the arrayer prints 16 spots and the print trays that where used
  # contains 384 wells. We guess that the printing was done during 
  # four days and each day 4 print trays where printed, i.e. 4*384=1536
  # spots. We also notice that the first day only 3.25 trays where printed,
  # i.e. 1248 spots. These findings has been verified with the lab that
  # produced the slides.
  days <- 1:4
  breakpoints <- 1248 + (days-1)*1536
  abline(v=breakpoints, col="gray")

  if (!exists("po.fitted")) {
    if (exists("ma.norm")) rm(ma.norm)
    ma.norm <- clone(ma)
  
  	# We decide to normalize for the print-order effects by normalizing each
  	# section (as specified by breakpoints) seperately. For each section we
  	# will estimate a smooth function and adjust the values according to
  	# this function. We will also scale the data in a way that all the
  	# sections has the same scale. This procedure is done both for the
  	# M and the A values.
  	po.fitted <- normalizePrintorder(ma.norm, breakpoints=breakpoints)
  }
  lines(as.vector(toPrintorderMatrix(layout, po.fitted$M[,1])))
  
  plotPrintorder(ma.norm, slide=1, breakpoints=breakpoints)

  Device$print(imgfile)
}
pushTag(out, "center")
printTag(out, "img", src=getName(imgfile))
popTag(out)
rm(imgfile)

#debug(normalizeSpatial.MicroarrayData)
#debug(normalizePrintorder.MicroarrayData)

########################################################################
# 2. Spatial artifacts
########################################################################
printTag(out, "h1", "Spatial artifacts")
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)
  subplots(4)

  plotSpatial(ma, slide=1)
  plot3d(ma, slide=1)

  # Next we want to normalize for spatial effects in both M and A.
  fitted <- normalizeSpatial(ma.norm)

  plotSpatial(ma.norm, slide=1)
  plot3d(ma.norm, slide=1)

  ma.tmp <- clone(ma.norm)
  smooth <- normalizeSpatial(ma.tmp, slide=1)
}

if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)

  xy <- getPosition(layout)
  x <- xy[,1]; x <- c(x,x)
  y <- xy[,2]; y <- c(y,y)
  z <- c(fitted$M[,1], smooth$M[,1])
  col <- as.vector(matrix(c("red", "blue"), nrow=nbrOfSpots(ma), 
                   ncol=2, byrow=TRUE))
  plot3d(x,-y,z, pch=".", col=col, phi=40, theta=30, 
         xlab="x", ylab="y", zlab="M")
}

########################################################################
# 2. Intensity dependent artifacts
########################################################################
printTag(out, "h1", "Intensity dependent artifacts")
if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)
  subplots(4)

  plot(ma, slide=1); lowessCurve(ma, gridwise=TRUE)
  boxplot(ma, slide=1)

  # Next we want to normalize for spatial effects in both M and A.
  normalizeWithinSlide(ma.norm, method="s")

  plot(ma.norm, slide=1); lowessCurve(ma.norm, gridwise=TRUE)
  boxplot(ma.norm, slide=1)
}



if (!Device$isOpen(fig <- fig + 1)) {
  Device$set(fig)
  subplots(4)

  if (!exists("bma"))
    bma <- as.BMAData(ma.norm)

  # The spots with the top 1% B values
  extreme <- top(bma$B, 0.01)

  down <- extreme & (bma$M < 0)
  up <- extreme & (bma$M > 0)

  plot(bma)
  highlight(bma, down, col="green"); highlight(bma, up, col="red")
  plot(bma, "BvsM")
  highlight(bma, down, col="green"); highlight(bma, up, col="red")
  plotSpatial(bma)
  highlight(bma, down, col="green"); highlight(bma, up, col="red")
  col <- getColors(bma); col[down] <- "green"; col[up] <- "red"
  pch <- rep(charToInt("."), nbrOfSpots(bma)); pch[extreme] <- 20
  plot3d(bma, col=col, pch=pch)
}


popTags(out)
printTimestamp(out)
close(out)

System$openBrowser(file)

