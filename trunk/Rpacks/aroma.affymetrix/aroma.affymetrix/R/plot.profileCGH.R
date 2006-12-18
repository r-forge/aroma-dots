setMethodS3("plot", "profileCGH", function(x, ylim=c(-1,1)*3, units="Mb", Bkp=FALSE, Smoothing="Smoothing", cnLevels=c(1/2,1,3/2), colDAGLAD=NULL, ticksBy=1, ...) {
  # To please R CMD check.
  this <- x;

  # Argument 'units':
  units <- match.arg(units);

  # Argument 'colDAGLAD':
  if (is.null(colDAGLAD)) { 
    colDAGLAD <- brewer.pal(5, "Dark2");
  }

  # Get the last physical position
  lastPos <- max(this$profileValues$PosBase);

  # X scale
  if (units == "Mb") {
    unit <- 6;
  }
  scale <- 10^unit;

  # Rescale
  lastPos <- lastPos / scale;
#  for (ff in c("profileValues", "profileValuesNA", "BkpInfo")) {
#    this[[ff]]$PosBase <- this[[ff]]$PosBase / scale;
#  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Plot GLAD fit
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  plotProfile(this, unit=unit, Bkp=Bkp, Smoothing=Smoothing,
                    colDAGLAD=colDAGLAD, ylim=ylim, ...);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add ruler
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (kk in 1:3) {
    at <- seq(from=0, to=lastPos, by=ticksBy*c(1,5,10)[kk]);
    tcl <- c(0.2,0.4,0.6)[kk];
    lwd <- c(1,1,2)[kk];
    for (ss in c(1,3))
      axis(side=ss, at=at, tcl=tcl, lwd=lwd, labels=FALSE);
  }
  text(x=at, y=par("usr")[3]-0.1, labels=at, srt=90, 
                                       adj=1, cex=1, xpd=TRUE);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add horizontal rulers
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  for (level in cnLevels) {
    abline(h=log2(level), col="blue", lty=2);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add number of data points
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  n <- nrow(this$profileValues);
  stext(text=sprintf("n=%d", n), side=4, pos=0, line=0, cex=0.7, col="lightgray");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Add std dev estimates
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  sdEst <- this$SigmaC$Value;
  text <- substitute(hat(sigma)==x, list=list(x=sprintf("%.3g", sdEst)));
  stext(text=text, side=3, pos=0.5, line=-2);
})




############################################################################
# HISTORY:
# 2006-11-29
# o Added std dev estimation to plot.
# 2006-10-30
# o Created from glad.R script from 2006-08-03.
############################################################################
