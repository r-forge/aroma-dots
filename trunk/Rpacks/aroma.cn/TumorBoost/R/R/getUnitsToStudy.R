getUnitsToStudy <- function(ces, chromosome, subset=NULL, ...) {
  # Get annotation data files
  cdf <- getCdf(ces);
  ugp <- AromaUgpFile$byChipType(getChipType(ces));

  unitTypes <- getUnitTypes(cdf);

  units <- getUnitsOnChromosome(ugp, chromosome);
  units <- units[unitTypes[units] == 2];

  # Subset?
  if (!is.null(subset)) {
    if (subset > 1) {
      n <- subset;
    } else {
      n <- subset*length(units);
    }
    idxs <- seq(from=1, to=length(units), length.out=n);
    idxs <- unique(as.integer(idxs));
    units <- units[idxs];
    rm(idxs);
  }

  pos <- getPositions(ugp, units=units);

  # Extract (thetaA,thetaB) for the tumor and normal
  thetaR <- extractTheta(ceNR, units=units, drop=TRUE);
  thetaR <- CartesianSnpData(thetaR);
  betaR <- asTotalFracBSnpData(thetaR);
  polarR <- asPolarSnpData(thetaR);

  # Extract (thetaA,thetaB) for the tumor and normal
  theta <- extractTheta(ces, units=units);

#  thetaN <- theta;
#  thetaN[,1,] <- normalizeAffine(theta[,1,], constraint=0.001)
#  thetaN[,2,] <- normalizeAffine(theta[,2,], constraint=0.001)
#  thetaN[,,1] <- thetaN[,,1] - 50;
#  thetaN[,1,] <- thetaN[,1,] + 1000;
#  str(thetaN);
#  theta <- thetaN;

  thetaT <- CartesianSnpData(theta[,,1]);
  thetaN <- CartesianSnpData(theta[,,2]);

  betaT <- asTotalFracBSnpData(thetaT);
  betaN <- asTotalFracBSnpData(thetaN);

  polarT <- asPolarSnpData(thetaT);
  polarN <- asPolarSnpData(thetaN);

  # Call normal genotypes
  betaNC <- callGenotypes(betaN);
  polarNC <- asPolarSnpData(betaNC);
  thetaNC <- asCartesianSnpData(betaNC);
  genotypesN <- as.integer(1+2*betaNC[,2]);

  # TumorBoost
  thetaTB <- pairedBoost(thetaT, thetaN);
  betaTB <- pairedBoost(betaT, betaN);
  polarTB <- pairedBoost(polarT, polarN);

  # Copy numbers
  C <- 2*betaT[,1]/betaN[,1];
  thetaC <- 2*asTotalFracBSnpData(thetaTB)[,1]/betaN[,1];
  betaC <- 2*betaTB[,1]/betaN[,1];
  polarC <- 2*asTotalFracBSnpData(polarTB)[,1]/betaN[,1];
  CN <- betaN[,1]/betaR[,1];

  data <- list(
    chromosome = chromosome,
    subset = subset,
    nbrOfUnits = length(units),
    unit = units, 
    position = pos/1e6,
    thetaR = thetaR,
    thetaT = thetaT,
    thetaTB = thetaTB,
    thetaN = thetaN,
    thetaNC = thetaNC,
    betaR = betaR,
    betaT = betaT,
    betaTB = betaTB,
    betaN = betaN,
    rho = betaT[,2]-betaN[,2],
    betaNC = betaNC,
    polarR = polarR,
    polarT = polarT,
    polarTB = polarTB,
    polarN = polarN,
    polarNC = polarNC,
    genotypesN = genotypesN,
    C = C,
    thetaC = thetaC,
    betaC = betaC,
    polarC = polarC,
    M = log2(C),
    CN = 2*CN,
    MN = log2(CN)
  );

  data;
} # getUnitsToStudy()

############################################################################
# HISTORY:
# 2009-03-31
# o Added getUnitsToStudy.R.
############################################################################

