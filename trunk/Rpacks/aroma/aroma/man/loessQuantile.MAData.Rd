%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  MAData.R
% 
% on Tue Jan 15 18:33:22 2008.
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{loessQuantile.MAData}
\alias{loessQuantile.MAData}
\alias{MAData.loessQuantile}
\alias{loessQuantile.MAData}
\alias{loessQuantile,MAData-method}


 \title{Estimates the (25\%,50\%,75\%) quantile M(A) curves using loess}

 \usage{\method{loessQuantile}{MAData}(this, slides=NULL, subset=NULL, method=c("each", "average"), factor=1.4826, family="symmetric", ...)}

 \arguments{
   \item{slides}{The slides which should be included in the calculations.
     If \code{\link[base]{NULL}}, all slides are included.}
   \item{subset}{Indices of spots to be used for the estimates. Selecting
     a random subset may speed up the estimation and still give the same
     result. If \code{\link[base]{NULL}}, all spots are included.}
   \item{method}{If \code{"each"}, estimates of the quantiles are given
     for each of the slides. If \code{"average"}, an estimate of the
     quantiles is given for all spots as if they came from the same slide.}
   \item{factor}{\code{\link[base]{numeric}} factor to multiply to the upper and lower
     quantiles. If the data is Gaussian distributed the default factor
     \code{1.4826} makes the upper and lower quantile to correspond to
     one standard deviation.}
   \item{family}{Family to be used by \code{\link[modreg]{loess}}.
     By default, a robust family is used.}
   \item{...}{Other arguments passed to \code{loess()}.}
 }

 \description{
   Estimates the (25\%,50\%,75\%) quantile M(A) curves using loess.
   This method is used by \code{\link[aroma:getAdjustedSpotVariability.MAData]{*getAdjustedSpotVariability}()}.
 }

 \examples{
   SMA$loadData("mouse.data")
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

   setMethodS3("lines", "LoessQuantile", function(object, ...) {
     lines(object$x, object$centre, ...);
     lines(object$x, object$above, ...);
     lines(object$x, object$below, ...);
   })

   # Four plots
   subplots(4)
   Alim <- Mlim <- NA
   for (k in 1:length(ma)) {
     Alim <- range(c(Alim, ma[[k]]$A), na.rm=TRUE)
     Mlim <- range(c(Mlim, ma[[k]]$M), na.rm=TRUE)
   }
 
   for (k in 1:length(ma)) {
     plot(ma[[k]], slide=1, xlim=Alim, ylim=Mlim, col=k+1)
     fit <- loessQuantile(ma[[k]], slide=1)[[1]]
     lines(fit, lwd=2)
   }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \seealso{
   \code{\link[aroma:getAdjustedSpotVariability.MAData]{*getAdjustedSpotVariability}()}.
   For more information see \code{\link{MAData}}.
 }
\keyword{internal}
\keyword{methods}
