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

