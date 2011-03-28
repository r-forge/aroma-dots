%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  backtransformAlleleRatiosByCentroids.matrix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 \name{backtransformAlleleRatiosByCentroids.matrix}
\alias{backtransformAlleleRatiosByCentroids.matrix}
\alias{backtransformAlleleRatiosByCentroids.matrix}

 \alias{backtransformAlleleRatiosByCentroids}

 \title{Reverse transform of angular signals}

 \description{
  Reverse transform of angular signals.
 }

 \usage{\method{backtransformAlleleRatiosByCentroids}{matrix}(beta, mu, truncate=FALSE, ...)}

 \arguments{
  \item{beta}{An JxI \code{\link[base]{numeric}} \code{\link[base]{matrix}} consisting of fractional signals for J SNPs across I samples.}
  \item{mu}{An JxC \code{\link[base]{numeric}} \code{\link[base]{matrix}} of fractional locations for C=2 or C=3 centroids.}
  \item{truncate}{If \code{\link[base:logical]{TRUE}}, the fractional signals are truncated to [0,1].}
  \item{...}{Not used.}
 }

 \value{
   Returns an JxI \code{\link[base]{numeric}} \code{\link[base]{matrix}} of backtransformed fractional signals.
 }

 \section{B Allele Fractions by Peiffer et al. (2007)}{
   In [2], the authors describe how to normalize B Allele Fractions (BAFs)
   based on predefined linear interpolation.  This method, which
   is used in the Illumina's BeadStudio software, is also described in [1].
   That method is obtained by passing \eqn{\theta} of [1] and [2] as
   argument \code{beta} while using a Jx3 \code{mu} centroid matrix and
   \code{truncate=TRUE}.
 }

 \examples{\dontrun{
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Example: Illumina HumanHap550 SNP data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Predefined genotype centroids for 5 SNPs
mu <- matrix(c(
  0.0042, 0.7057, 0.9933,
  0.0001, 0.5039, 0.9929,
  0.0195, 0.7133, 0.9652,
  0.0111, 0.4507, 0.9850,
  0.0036, 0.5049, 0.9817
), nrow=5, ncol=3, byrow=TRUE);
colnames(mu) <- c("AA", "AB", "BB");
rownames(mu) <- c("rs10013734", "rs10013737", "rs10013782", "rs10013797", "rs10013819");


# Some data
theta <- matrix(c(
  -0.0042,
  0.7057,
  0.9933,
  0.5039,
  0.0111
), nrow=5, ncol=1, byrow=TRUE);

# B Allele Fractions (BAFs) according to Peiffer et al. (2007)
thetaN <- backtransformAlleleRatiosByCentroids(theta, mu=mu, truncate=TRUE)


 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] Illumina,
      \emph{BeadStudio Genotyping Module v3.2 - User Guide},
      Illumina Inc., Part \#11284301, Rev. A, 2007, pp 114--115.\cr
  [2] Peiffer et al.,
      \emph{High-resolution genomic profiling of chromosomal aberrations
      using Infinium whole-genome genotyping},
      Genome Research, 2006, 16.\cr
 }
\keyword{methods}
