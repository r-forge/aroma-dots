#########################################################################/**
# @RdocClass TMAData
#
# @title "The TMAData class"
#
# \description{
#  \bold{IMPORTANT: This class is in a beta version and it might be 
#  replaced by something else.}\cr
#
#  @classhierarchy
#
#  Creates a new \code{TMAData} object.
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab \code{M} \tab The average M, i.e. the log of ratios between R and G,
#      over all slides. \cr
#    \tab \code{A} \tab The average A, i.e. the log of the intensities 
#      of R and G, over all slides. \cr
#    \tab \code{T} \tab The T values for each spot. \cr
#    \tab \code{df} \tab The \emph{degrees of freedom} for each spot. \cr
#    \tab \code{SE} \tab The standard error for each spot. \cr
#  }
#
#  @allmethods
# }
#
# \details{
#   Given a set of slides stored in  a \link{MAData} object, which has
#   been normalized, the "average" slide can be calculated. The average
#   slide \eqn{(M',A')} is calculated (spot wise) as:
#    \deqn{M'=\sum_{i=1}{N}M_i/N}{M'=(1/N)*sum of of a M}, and
#    \deqn{A'=\sum_{i=1}{N}A_i/N}{A'=(1/N)*sum of of a A}.
#
#   In addition to M' and A', the T values are calculated. The T value
#   for a spot is the average M for that spot across all slides divided by
#   their standard error. The degrees of freedom is the number of slides
#   included when calculated the T values. For most spots, this is the
#   same as the number slides, but if some spots is registered as 
#   @NA they will be excluded from the calculation of the T value
#   and the degrees of freedom will be smaller.
# }
#
# \note{
#   Note that the standard error, SE, is a redundant field. It can easily be
#   calculated as \eqn{SE=M/T} or more explicite as \code{SE=tma$M/tma$T}.
#   However, since it might be used very often it is included.
# }
#
# @author
#
# \examples{
#   SMA$loadData("mouse.data")
#   layout <- Layout$read("MouseArray.Layout.dat", path=system.file("data-ex", package="aroma"))
#   raw <- RawData(mouse.data, layout=layout)
#   ma <- getSignal(raw, bgSubtract=TRUE)
#   normalizeWithinSlide(ma, "s")
#   normalizeAcrossSlides(ma)
#   tma <- as.TMAData(ma)
#
#   layout(matrix(1:4, ncol=2, byrow=TRUE))
#   plot(tma, "MvsA")
#   plot(tma, "TvsA")
#   plot(tma, "TvsSE")
#   plot(tma, "spatialM")
# }
#
# \seealso{
#   @see "MAData.as.TMAData".
# }
#*/#########################################################################
setConstructorS3("TMAData", function(M=NULL, A=NULL, T=NULL, df=NULL, stderr=NULL, layout=NULL, extras=list()) {
  # SE is redundant data, but simplifies life. The memory is not an issue either,
  # at least not if one compares to the huge GPR objects or large MAData objects.
  if (!is.null(M))           M <- as.matrix(M);
  if (!is.null(A))           A <- as.matrix(A);
  if (!is.null(T))           T <- as.matrix(T);
  if (!is.null(df))         df <- as.matrix(df);
  if (!is.null(stderr)) stderr <- as.matrix(stderr);

  this <- extend(MAData(M=M, A=A, layout=layout, extras=extras), "TMAData",
    T      = T,
    df     = df,
    stderr = stderr
  )

  this$.fieldNames <- c("M", "A", "T", "df", "SE");
  
  this <- setLabel(this, "T",  expression(T==bar(M)/SE(M)));
  this <- setLabel(this, "SE", expression(SE(M)));
  this <- setLabel(this, "M",  expression(bar(M)));
  this <- setLabel(this, "A",  expression(bar(A)));
  this;
}, deprecated=TRUE)



############################################################################
############################################################################
## 
##  DATA STRUCTURAL METHODS
## 
############################################################################
############################################################################


setMethodS3("as.TMAData", "TMAData", function(this, ...) {
  this;
}, deprecated=TRUE)


setMethodS3("as.TMAData", "ANY", function(this, ...) {
  TMAData(this, ...)
}, deprecated=TRUE)


setMethodS3("as.character", "TMAData", function(this) {
  s <- data.class(this);
  s <- paste(sep="",s,": T ",      com.braju.sma.dimStr(this$T));
  s <- paste(sep="",s,", df ",     com.braju.sma.dimStr(this$df));
  s <- paste(sep="",s,", stderr ", com.braju.sma.dimStr(this$stderr));
  s <- paste(sep="",s,", ",        as.character.MAData(this));
  s;
}, deprecated=TRUE)


############################################################################
############################################################################
## 
##  PLOTTING & GRAPHICAL METHODS
## 
############################################################################
############################################################################
# Generate grayscale colors by default.
setMethodS3("getColors", "TMAData", function(this, what="T", slide=1, log=NULL, palette=NULL, ...) {
#  if (!is.null(what) && !is.element(what, getFieldNames(this)))
#    throw("Argument 'what' is refering to an unknown field: ", what);

  # If M vs A is requested use getColors() of MAData...
  if ("M" %in% what && "A" %in% what)
    return(getColors(as.MAData(this)));

  if (any(what == "T"))
    what <- what[what == "T"]
  else
    what <- what[1];
  x <- this[[what]][,slide];
  if (!is.null(log))
    x <- log(x, base=log);
  dim.range <- c(0,0.8);
  if (is.null(palette) || palette == "auto")
    palette <- "heat";
  if (palette == "grayscale")
    Colors$getGray(x, dim.range=dim.range)
  else if (palette == "heat")
    Colors$getHeatColors(x, dim.range=dim.range);
}, deprecated=TRUE)

setMethodS3("plotXY", "TMAData", function(this, what=c("SE","T"), ...) {
  plotXY.MAData(this, what=what, ...);
}, deprecated=TRUE)


setMethodS3("plot", "TMAData", function(x, what="TvsA", ...) {
  # To please R CMD check...
  this <- x;

  plot.MAData(this, what=what, ...);
}, deprecated=TRUE)

setMethodS3("plotSpatial", "TMAData", function(this, what="T", ...) {
  plotSpatial.MAData(this, what=what, ...);
}, deprecated=TRUE)


setMethodS3("plotPrintorder", "TMAData", function(this, what="T", ...) {
  plotPrintorder.MAData(this, what=what, ...);
}, deprecated=TRUE)




############################################################################
############################################################################
## 
##  STATISTICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("topSpots", "TMAData", function(this, X="A", Y="T", limits=NULL, ...) {
  if (X == "pval" || Y == "pval")
    this$pval <- as.matrix(tma$getPValues()$pval);
  if (is.null(limits))
    if (Y == "pval") limits <- "lower" else limits <- "both";

  topSpots.MAData(this, X=X, Y=Y, limits=limits, ...);
}, protected=TRUE, deprecated=TRUE);


setMethodS3("getPValues", "TMAData", function(this, 
           alternative=c("two.sided", "less", "greater"), conf.level=0.95) {

  # Define T as:
  #
  #   T = mean(X)/SE(X) = mean(X)/(var(X)/n)
  #
  # The t distribution with df = n degrees of freedom has density 
  #
  #    f(x) = Gamma((n+1)/2) / (sqrt(n pi) Gamma(n/2)) (1 + x^2/n)^-((n+1)/2)
  #
  # for all real x. It has mean 0 (for n > 1) and variance n/(n-2)
  # (for n > 2). 
  #
  # Assuming that X (here M) is N(0, sigma^2), then T (here 'tstat'),
  # is distributed as non-centrally t with df= n-1 degrees of
  # freedom and non-centrality parameter ncp=0 (=0-0).

  T <- this$T[,1];
  df <- this$df[,1];

  if (alternative == "less") {
    pval <- pt(T, df);
    cint <- c(NA, T + qt(conf.level, df));
  } else if (alternative == "greater") {
    pval <- pt(T, df, lower = FALSE);
    cint <- c(T - qt(conf.level, df), NA);
  } else if (alternative == "two.sided") {
    pval <- 2 * pt(-abs(T), df);
    alpha <- 1 - conf.level;
    warn.def <- getOption("warn");
    options(warn=-1);               # Ignore warnings for 'qt' for
                                    # cases where df==0.
    cint <- qt(1 - alpha/2, df);
    options(warn=warn.def);
    cint <- T + c(-cint, cint);
  } else
    throw("Unknown alternative: ", alternative);
 
  list(pval=pval, cint=cint);
}, deprecated=TRUE)


setMethodS3("getDegreesOfFreedom", "TMAData", function(this) {
  this$df;
}, deprecated=TRUE)

setMethodS3("getSE", "TMAData", function(this) {
  cbind(this$stderr, this$stderr);
}, deprecated=TRUE)

setMethodS3("getStandardError", "TMAData", function(this) {
  getSE(this);
}, deprecated=TRUE)


############################################################################
# HISTORY:
# 2003-04-08
# o getColors() of TMAData is (once again) returning red-green colors if
#   colors for (M,A) is requested.
# 2002-10-09
# o Added getColors(), which now by default for T generate heat colors.
# o Update plotXY(), added plot(), plotSpatial(), plotPrintorder().
# o Removed redundant field SE and replaced it with getSE().
# 2002-05-11
# * BUG FIX: Forgot to assure that M and A are stored as matrices.
# 2002-02-26
# * Remove as.data.frame(), get() etc which are now generic in the class
#   MicroarrayData.
# * Updated the code to make use of setMethodS3().
# 2001-08-08
# * BUG FIX: Removed an erronous as.MAData.MAData that shouldn't be in this
#   file.
# 2001-08-01
# * Added the field spot in as.data.frame().
# 2001-07-15
# * Bug fix: Forgot 'eps' in equals(), which is now also comparing df & SE.
# 2001-07-14
# * Added the *redundant* field SE, which can be calculated as SE=M/T. 
#   However, I believe that it will be used so much so it should be included
#   by default. It would be pretty cool if one could create virtual fields
#   that are calculated on demand. Should be possible to detect and do
#   through "$<-.Object" etc. A future work though!
# * Update equals to include Layout comparison too. Added Rdoc for equals().
# 2001-07-12
# * Now plotting SE is also supported. Creates a temporary SE if needed.
# 2001-07-05
# * Updated topSpots() and plot.TMAData(). Most of the plot functionalities
#   are now in plot.MicroarrayData().
# 2001-07-04
# * Created getPValues().
# * Created the protected method plotYvsX() and the generic plot() method.
# 2001-04-05
# * Forgot to declare highlight = METHOD.
# 2001-04-04
# * Added TvsSE support for highlight.
# * Added set/getTLabel().
# 2001-04-02
# * plotTvsA(): Added default values if styler is not given.
# 2001-03-28
# * Created.
############################################################################
