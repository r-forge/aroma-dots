setConstructorS3("MicroarrayWeights", function(...) {
  extend(Object(), "MicroarrayWeights",
    ...
  )
})


#########################################################################/**
# @set "class=MicroarrayWeights"
# @RdocMethod "unite"
#
# @title "Coalesce weights of lower order into higher order weights"
#
# @synopsis
#
# \description{
#  @get "title". 
#
#  Each new weight is the square root of the average sum of squares of
#  lower-order weights.
# }
#
# \arguments{
#  \item{weights}{A NxK @matrix of weights to be united row by row.}
# }
#
# \value{
#   Returns a weight @vector of length N.
# }
#
# @author
#
# \seealso{
#   @seeclass.
# } 
#*/######################################################################### 
setMethodS3("unite", "MicroarrayWeights", function(static, weights) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'weights'
  if (!is.matrix(weights))
    throw("Argument 'weights' must be a matrix.");

  K <- ncol(weights);
  weights <- sqrt(rowSums(weights^2) / K);

  weights;
}, static=TRUE)



setConstructorS3("ProbeWeights", function(...) {
  extend(MicroarrayWeights(...), "SignalWeights");
})


setConstructorS3("SignalWeights", function(...) {
  extend(MicroarrayWeights(...), "SignalWeights");
})

   

#########################################################################/**
# @set "class=SignalWeights"
# @RdocMethod "fromSaturatedSignals"
#
# @title "Weight function for scanner signals"
#
# @synopsis
#
# \description{
#  @get "title". 
# }
#
# \arguments{
#  \item{x}{A @vector (or a @matrix) of scanner signals for which weights
#    are calculated.}
#  \item{satRatio}{A @double in [0,1] specifying the relative level (of the 
#    maximum signal) where saturations starts.}
#  \item{maxSignal}{A @numeric specifying the maximum scanner signal.}
# }
#
# \value{
#   Returns a @vector (or a @matrix) of weights in [0,1].
# }
#
# \details{
#  All signals less than \code{satRatio*maxSignal} get weight one, all 
#  greater or equal to \code{maxSignal} get weight zero, and inbetween
#  there is a linear decrease from one to zero.
#
#  Default value of argument \code{satRatio} is chosen such that saturation
#  "starts" about 50000 for signals on the range [0,65535], which was
#  reported by [1].
# }
#
# \section{Missing values}{
#   Missing-value signals (@NA), get zero weight.
# }
#
# \examples{
#   x <- seq(from=0, to=10^5, by=100)
#   plot(x, SignalWeights$fromSaturatedSignals(x), type="l")
# }
#
# @author
#
# \references{
#  [1] H. Lyng, A. Badiee, D.H. Svendsrud, E. Hovig, O. Myklebost, T. Stokke. 
#      Profound influence of microarray scanner characteristics on gene 
#      expression ratios: analysis and procedure for correction. 
#      BMC Genomics, 2004, 5:10.
# }
#
# \seealso{
#   @seeclass.
# } 
#*/######################################################################### 
setMethodS3("fromSaturatedSignals", "SignalWeights", function(static, x, satRatio=0.75, maxSignal=2^16-1) { 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'x':
  if (any(x < 0, na.rm=TRUE))
    throw("Argument 'x' must be non-negative (or NA).");

  # Argument 'satRatio':
  if (length(satRatio) != 1 || !is.numeric(satRatio))
    throw("Argument 'satRatio' must be a scalar in [0,1].");
  if (satRatio < 0 || satRatio > 1)
    throw("Argument 'satRatio' is out of range [0,1]: ", satRatio);

  # Argument 'maxSignal':
  if (length(maxSignal) != 1 || !is.numeric(maxSignal))
    throw("Argument 'maxSignal' must be a positive scalar.");
  if (maxSignal <= 0)
    throw("Argument 'maxSignal' must be positive: ", maxSignal);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Weight function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  w <- 0*x;

  # Missing signals get zero weight
  nas <- is.na(x);
  w[nas] <- 0;

  nonSaturated <- (!nas & x < satRatio*maxSignal);
  w[nonSaturated] <- 1;
  saturated <- (!nas & !nonSaturated & x < maxSignal);
  w[saturated] <- (maxSignal-x[saturated])/((1-satRatio)*maxSignal); 


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Return weights
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  w;
}, static=TRUE)


############################################################################
# HISTORY:
# 2006-02-08
# o Rd examples used non-existing function.
# o Rd bug fix: Replaced section 'example' with 'examples'.
# 2005-02-08
# o Created. Sunny notes written yesterday at Bar Italia in Västra Hamnen.
############################################################################
