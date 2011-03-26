###########################################################################/**
# @set "class=matrix"
# @RdocMethod backtransformAlleleRatiosByCentroids
# @alias backtransformAlleleRatiosByCentroids
#
# @title "Reverse transform of angular signals"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#  \item{beta}{An JxI @numeric @matrix consisting of fractional signals for J SNPs across I samples.}
#  \item{mu}{An JxC @numeric @matrix of fractional locations for C=2 or C=3 centroids.}
#  \item{truncate}{If @TRUE, the fractional signals are truncated to [0,1].}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns an JxI @numeric @matrix of backtransformed fractional signals.
# }
#
# \section{B Allele Fractions by Peiffer et al. (2007)}{
#   In [2], the authors describe how to normalize B Allele Fractions (BAFs)
#   based on predefined piece-wise linear regressions.  This method, which
#   is used in the Illumina's BeadStudio software, is also described in [1].
#   That method is obtained by passing \eqn{\theta} of [1] and [2] as
#   argument \code{beta} while using a Jx3 \code{mu} centroid matrix and
#   \code{truncate=TRUE}.
# }
#
# \examples{\dontrun{
#  @include "../incl/backtransformAlleleRatiosByCentroids.Rex"
# }}
#
# @author
#
# \references{
#  [1] Illumina, 
#      \emph{BeadStudio Genotyping Module v3.2 - User Guide}, 
#      Illumina Inc., Part \#11284301, Rev. A, 2007, pp 114--115.\cr
#  [2] Peiffer et al., 
#      \emph{High-resolution genomic profiling of chromosomal aberrations
#      using Infinium whole-genome genotyping}, 
#      Genome Research, 2006, 16.\cr
# }
#*/###########################################################################
setMethodS3("backtransformAlleleRatiosByCentroids", "matrix", function(beta, mu,  truncate=FALSE, ...) {
  # Argument 'beta':
  dim <- dim(beta);

  # Argument 'mu':
  stopifnot(is.matrix(mu));
  stopifnot(nrow(mu) == dim[1]);
  nbrOfCentroids <- ncol(mu);
  stopifnot(nbrOfCentroids == 2 || nbrOfCentroids == 3);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (a) Adjust such that AA and BB centroids are at 0 and 1
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  a <- mu[,1];
  b <- mu[,nbrOfCentroids] - mu[,1];

  beta <- beta - a;
  beta <- beta / b;

  mu <- mu - a;
  mu <- mu / b;

  # Sanity check
  stopifnot(all(mu[,1] == 0));
  stopifnot(all(mu[,nbrOfCentroids] == 1));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (b) Adjust such that AB centroids are at 1/2
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (nbrOfCentroids == 3) {
    c <- mu[,2] / (1/2);
    d <- (1-mu[,2]) / (1/2);
    for (kk in seq(length=ncol(beta))) {
      keep <- (beta[,kk] < mu[,2]);
      keep <- keep & is.finite(keep);
      beta[keep,kk] <- 0 + (beta[keep,kk] - 0) / c[keep];
      keep <- !keep;
      beta[keep,kk] <- 1 - (1 - beta[keep,kk]) / d[keep];
    }

    # Same for centroids; just for validation
    c <- mu[,2] / (1/2);
    d <- (1-mu[,2]) / (1/2);
    for (kk in seq(length=nbrOfCentroids)) {
      keep <- (mu[,kk] < mu[,2]);
      keep <- keep & is.finite(keep);
      mu[keep,kk] <- mu[keep,kk] / c[keep];
      keep <- !keep;
      mu[keep,kk] <- 1 - (1 - mu[keep,kk]) / d[keep];
    }

    # Sanity check
    stopifnot(all(mu[,2] == 1/2));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # (c) Truncate beta to [0,1]
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (truncate) {
    beta[beta < 0] <- 0;
    beta[beta > 1] <- 1;
  }

  attr(beta, "muN") <- mu;

  beta;
}) # backtransformAlleleRatiosByCentroids()



############################################################################
# HISTORY:
# 2011-03-25
# o Added backtransformAlleleRatiosByCentroids().
# o Created.
############################################################################
