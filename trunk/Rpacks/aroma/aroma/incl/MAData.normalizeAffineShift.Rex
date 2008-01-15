SMA$loadData("mouse.data")
layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
raw <- RawData(mouse.data, layout=layout)
methods <- c("Kerr", "Bengtsson")
Alim <- c(6,16)
Mlim <- c(-1,1)*3
subplots(2*length(methods))
for (method in methods) {
  ma <- getSignal(raw, bgSubtract=TRUE)
  fit <- normalizeAffineShift(ma, method=method, slide=1)
  a <- fit$shift[1,1]
  msg <- sprintf("Optimal shift according to '", method, "' is a = %1.f.\n", a)
  cat(msg)
  plot(ma, xlim=Alim, ylim=Mlim)
  lowessCurve(ma)
  
  x <- seq(-2*a,2*a, length=200)
  y <- sapply(x, FUN=fit$objectiveFunction, R=getR(ma),G=getG(ma))
  plot(x,y, pch=176, xlab="shift", ylab="objective");
  lines(x,y)
  abline(v=a)
} # for (method ...)
