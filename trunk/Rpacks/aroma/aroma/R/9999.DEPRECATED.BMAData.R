#########################################################################/**
# @RdocClass BMAData
#
# @title "The BMAData class"
#
# \description{
#  \bold{IMPORTANT: This class is in a beta version and it might be 
#  replaced by something else.}\cr
#
#  @classhierarchy
#
#  Creates a new \code{BMAData} object.
#
#  \bold{Note:} We recommend that you use Gordon Smyth's limma package 
#  (\url{http://bioinf.wehi.edu.au/limma/}) and the @see "limma::ebayes" 
#  function instead, which is also based on [1]. Moreover, the limma package 
#  is richly documented. /July, 2003.
# }
#
# \section{Fields and Methods}{
#  \bold{Fields}
#  \tabular{rll}{
#    \tab \code{M} \tab The average M, i.e. the log of ratios between R and G,
#      over all slides. \cr
#    \tab \code{A} \tab The average A, i.e. the log of the intensities 
#      of R and G, over all slides. \cr
#    \tab \code{B} \tab The B values for each spot. \cr
#  }
#
#  @allmethods
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
#   bma <- as.BMAData(ma, p=0.01)
#
#   fM <- MFilter(bma, top=0.01, col="red")
#   fB <- BFilter(bma, top=0.01, col="black")
#   fAnd <- AndFilter(fM, fB, col="green")
#
#   layout(matrix(1:4, ncol=2, byrow=TRUE))
#   plot(bma, "MvsA"); highlight(fAnd, recursive=TRUE);
#   plot(bma, "BvsM"); highlight(fAnd, recursive=TRUE);
#   plot(bma, "BvsA"); highlight(fAnd, recursive=TRUE);
# }
#
# \references{
#  [1] Lönnstedt, I. and Speed, T. P. (2002). Replicated microarray data. 
#      \emph{Statistica Sinica} \bold{12}, 31-46. 
# }
#
# \seealso{
#   @see "limma::ebayes" (\url{http://bioinf.wehi.edu.au/limma/}).
#   @see "MAData.as.BMAData". Compare the results with the
#   \emph{t}-statistics method @see "TMAData".
# }
#*/#########################################################################
setConstructorS3("BMAData", function(M=NULL, A=NULL, B=NULL, df=NULL, params=NULL, layout=NULL, extras=list()) {
  # SE is redundant data, but simplifies life. The memory is not an issue either,
  # at least not if one compares to the huge GPR objects or large MAData objects.
  if (!is.null(M)) M <- as.matrix(M);
  if (!is.null(A)) A <- as.matrix(A);
  if (!is.null(B)) B <- as.matrix(B);

  this <- extend(MAData(M=M, A=A, layout=layout, extras=extras), "BMAData", 
    B      = B,
    df     = df,
    params = params
  )

  this$.fieldNames <- c("M", "A", "B", "df");
  
  this <- setLabel(this, "M",  expression(bar(M)));
  this <- setLabel(this, "A",  expression(bar(A)));
  this <- setLabel(this, "B",  expression(B));
  this;
}, deprecated=TRUE)


############################################################################
############################################################################
## 
##  DATA STRUCTURAL METHODS
## 
############################################################################
############################################################################


setMethodS3("as.BMAData", "BMAData", function(this, ...) {
  this;
}, deprecated=TRUE)

setMethodS3("as.BMAData", "ANY", function(this, ...) {
  BMAData(this, ...);
}, deprecated=TRUE)

setMethodS3("as.character", "BMAData", function(this) {
  s <- data.class(this);
  s <- paste(sep="",s,": B ", com.braju.sma.dimStr(this$B));
  s <- paste(sep="",s,", df ", this$df);
  s <- paste(sep="",s,", parameters (", paste(names(this$params),
                         this$params, sep="=", collapse=", "), ")");
  s <- paste(sep="",s,", ", as.character.MAData(this));
  s;
}, deprecated=TRUE)



############################################################################
############################################################################
## 
##  PLOTTING & GRAPHICAL METHODS
## 
############################################################################
############################################################################

setMethodS3("plotXY", "BMAData", function(this, what=c("A","M"), ...) {
  plotXY.MAData(this, what=what, ...);
}, deprecated=TRUE)



############################################################################
############################################################################
## 
##  STATISTICAL METHODS
## 
############################################################################
############################################################################



############################################################################
# HISTORY:
# 2003-07-07
# o Added links and recommendations about the limma package instead.
# 2002-02-26
# * Removed as.data.frame() and get(), which are now generic in the
#   MicroarrayData class.
# Update to make use of setMethodS3().           
# 2001-09-29
# * Created from TMAData.
############################################################################


